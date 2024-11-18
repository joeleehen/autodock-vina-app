#!/usr/bin/env python3

# NOTE: I added changes to final_sort for benchmarking purposes
# THIS FILE HAS CHANGES THAT SHOULD NOT BE MERGED TO MAIN!

import argparse
import logging
import os
from os.path import exists
from os.path import basename
import pickle # For unpickling pickled ligand files
import shutil
import subprocess
import time

import blosc # For decompressing compressed ligand files
from mpi4py import MPI
from vina import Vina


# Setup base MPI declarations
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

# Parse user inputs
parser = argparse.ArgumentParser(description ='vina commands')
parser.add_argument('-r', '--receptor',type=str,
                    help='Receptor for processing')
parser.add_argument('-c', '--center',type=str, required=True,
                    help='center position')
parser.add_argument('-s', '--size', type=str, required=True,
                    help='box size')
parser.add_argument('-f', '--forcefield', type=str, required=True,
                    help='vina or ad4 forcefield')
parser.add_argument('-d', '--docking', type=str, required=True,
                    help='basic or flexible docking')
parser.add_argument('-x', '--sidechains',type=str, default='EMPTY', required=False,
                    help='sidechains if flexible docking')
parser.add_argument('-n', '--number', type=int, required=True,
                    help='top n results to return')
parser.add_argument('-l', '--ligand_library', type=str, required=True,
                    help='ligand library')
args = parser.parse_args()

# Initialize logging
format_str=f'[%(asctime)s {RANK}] %(filename)s:%(funcName)s:%(lineno)s - %(levelname)s: %(message)s'
logging.basicConfig(level=logging.DEBUG, format=format_str)
#logging.basicConfig(level=logging.DEBUG, format=format_str, filename='autodock.log', filemode='w')

# Global constants
CENTER_X = float((args.center).split(',')[0])
CENTER_Y = float((args.center).split(',')[1])
CENTER_Z = float((args.center).split(',')[2])
SIZE_X = float((args.size).split(',')[0])
SIZE_Y = float((args.size).split(',')[1])
SIZE_Z = float((args.size).split(',')[2])
USER_CONFIGS = {'center_x': CENTER_X,
                'center_y': CENTER_Y,
                'center_z': CENTER_Z,
                'size_x': SIZE_X,
                'size_y': SIZE_Y,
                'size_z': SIZE_Z}
FULL_RECEPTOR=args.receptor
RECEPTOR=FULL_RECEPTOR.split('.')[0]
FLEX_RECEPTOR=f'{RECEPTOR}_flex'
DOCKING = args.docking
if DOCKING == 'rigid':
    FLEXIBLE = False
elif DOCKING == 'flexible':
    FLEXIBLE = True
SIDECHAINS = (args.sidechains).split('_')
LIBRARY_NAME = args.ligand_library.split('/')[4]
NUMBER_OF_OUTPUTS = args.number if args.number <= 1000 else 1000

logging.debug(f'CENTER = {args.center}; BOX = {args.size}; DOCKING = {DOCKING}')
logging.debug(f'LIGAND LIBRARY = {args.ligand_library}; LIBRARY NAME = {LIBRARY_NAME}')

# Internal constants
# tasks should be nodes * 128 / cpus
# These values were determined through internal benchmarks to allow an entire set to run
# under 24 hours
TASKS = int(os.environ['SLURM_NTASKS']) # What the user chose on the web portal
NODES = int(os.environ['SLURM_NNODES']) # What the user chose on the web portal
if LIBRARY_NAME in ['Test-set-compressed', 'Enamine-PC-compressed', 'ZINC-fragments-compressed', 'ZINC-in-trials-compressed']:
    EXPECTED_NODES = 1
    EXPECTED_TASKS = 32
elif LIBRARY_NAME == 'Enamine-HTSC-compressed':
    EXPECTED_NODES = 10
    EXPECTED_TASKS = 320
elif LIBRARY_NAME == 'Enamine-AC-compressed':
    EXPECTED_NODES = 3
    EXPECTED_TASKS = 96
CPUS = 4
VERBOSITY = 0 # Prints vina docking progress to stdout if set to 1 (normal) or 2 (verbose)
POSES = 1 # If set to 1, only saves the best pose/score to the output ligand .pdbqt file
EXHAUSTIVENESS = 8
MAX_SIDECHAINS = 6
MAX_BOXSIZE = 30
FILES_BEFORE_CLEANUP = 1000

logging.debug(f'TASKS = {TASKS}; NODES = {NODES}; EXPECTED TASKS = {EXPECTED_TASKS}; EXPECTED NODES = {EXPECTED_NODES}')



def pre_processing():
    '''
    This function is called by Rank 0 only. It prepares the config files, the 
    receptor, and the affinity maps. It also verifies the user configuration
    '''
    
    os.makedirs('configs')
    os.makedirs('output/pdbqt')
    os.makedirs('output/results/ligands')
    prep_config()
    prep_receptor()
    prep_maps()
    check_user_configs()

    return


def prep_config():
    '''
    Write user input configurations (grid center and box size) to a config file
    for AutoDock Vina to use
    '''

    with open('./configs/config.config', 'w+') as fout:
        for config, value in USER_CONFIGS.items():
            fout.write(f'{config} = {value}\n')
    logging.debug('Config file written')

    return


def prep_receptor():
    '''
    Converts a PDB receptor to a PDBQT receptor, if needed. If the user has
    specified flexible docking, also prepares the rigid receptor and user-chosen
    flexible sidechains. The prepare_receptor utility is part of the ADFR suite.
    '''

    if FULL_RECEPTOR.endswith('.pdb'):
        try:
            subprocess.run([f'prepare_receptor -r {RECEPTOR}.pdb'], shell=True)
            logging.debug('Receptor prepped successfully')
        except Exception as e:
            logging.error(f'Error prepping receptor')
            logging.error(e)
            COMM.Abort()

    if FLEXIBLE == True:
        try:
            subprocess.run([f"pythonsh /autodock-src/scripts/prepare_flexreceptor.py \
                                       -g {RECEPTOR}.pdbqt -r {RECEPTOR}.pdbqt \
                                       -s {'_'.join(SIDECHAINS)}"], shell=True)
            logging.debug('Flex receptor prepped successfully')
        except Exception as e:
            logging.error(f'Error prepping flex receptor')
            logging.error(e)
            COMM.Abort()

    return


def prep_maps():
    '''
    Generate affinity maps using the write-gpf.py utility script (which is part
    of this app) and ADFR's autogrid4. This process wil generate maps for all
    possible ligand atom types.
    '''

    if args.forcefield == 'ad4':
        if exists(f'{RECEPTOR}.gpf'):
            os.remove(f'{RECEPTOR}.gpf')
        subprocess.run([f"python3 /autodock-src/scripts/write-gpf.py --box {'./configs/config.config'} \
                        {RECEPTOR}.pdbqt"], shell=True)
        subprocess.run([f'autogrid4 -p {RECEPTOR}.gpf'], shell=True)
        logging.debug('Affinity maps generated')

    return


def check_user_configs():
    '''
    Checking for:
        Box size does not exceed MAX_BOXSIZE in any dimension
        Receptor is pdb or pdbqt file
        Grid center must be within receptor min/max bounds
        Sidechains must be in the pdbqt file
        Sidechains must not exceed MAX_SIDECHAINS (for flexible)
        Number of nodes and tasks must match expected

    Parsing coordinates according to:
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    '''
    # User input box size must be within bounds specified below
    for size in [SIZE_X, SIZE_Y, SIZE_Z]:
        if not (size <= MAX_BOXSIZE and size >= 1):
            logging.error(f'box size is outside the bounds (1-{MAX_BOXSIZE}).')
            COMM.Abort()

    # User must input a file ending in .pdb or .pdbqt
    if not (FULL_RECEPTOR.endswith('.pdb') or FULL_RECEPTOR.endswith('.pdbqt')):
        logging.error('Please provide a .pdb or .pdbqt file.')
        COMM.Abort()

    # Find bounds of pdbqt file in x,y,z dimensions
    all_sidechains = []
    xbounds = []
    ybounds = []
    zbounds = []
    with open(f'{RECEPTOR}.pdbqt', 'r') as fin:
        line = fin.readline()
        while line:
            line = fin.readline()
            if line.startswith('ATOM') or line.startswith('HETATM'):
                xbounds.append(float(line[30:38]))
                ybounds.append(float(line[38:46]))
                zbounds.append(float(line[46:54]))
                all_sidechains.append(line[17:20] + line[22:26])

    logging.debug(f'X-bounds = {min(xbounds)} - {max(xbounds)} | '
                  f'Y-bounds = {min(ybounds)} - {max(ybounds)} | '
                  f'Z-bounds = {min(zbounds)} - {max(zbounds)}')       

    # Confirm center is within bounds of protein
    if not (min(xbounds)) <= CENTER_X <= (max(xbounds)):
        logging.error(f'Center x coordinate {CENTER_X} is not within bounds.')
        COMM.Abort()
    if not (min(ybounds)) <= CENTER_Y <= (max(ybounds)):
        logging.error(f'Center y coordinate {CENTER_Y} is not within bounds.')
        COMM.Abort()
    if not (min(zbounds)) <= CENTER_Z <= (max(zbounds)):
        logging.error(f'Center z coordinate {CENTER_X} is not within bounds.')
        COMM.Abort()

    # Confirm sidechains exist
    if FLEXIBLE == True:
        for sidechain in SIDECHAINS:
            if not sidechain in all_sidechains:
                logging.error(f'Please provide valid flexible sidechain names '
                              f'separated by underscores (e.g. THR315_GLU268).')
                COMM.Abort()
    
    # Maximum of 6 sidechains
    if len(SIDECHAINS) > MAX_SIDECHAINS:
        logging.error(f'Too many sidechains specified (max: {MAX_SIDECHAINS}).')
        COMM.Abort()
    
    # User inputted #Nodes and #Tasks must match our internal values (specified above) exactly
    if not (TASKS == EXPECTED_TASKS) or not (NODES == EXPECTED_NODES):
        logging.error(f'Incorrect values for #Nodes and/or #ProcessorsPerNode.\n'
                      f'Please review input guidelines before submitting a job.\n'
                      f'Current #Nodes={NODES}.\n'
                      f'Current#Tasks={TASKS}.\n'
                      f'Expected #Nodes for {LIBRARY_NAME}={EXPECTED_NODES}.\n'
                      f'Expected #Tasks for {LIBRARY_NAME}={EXPECTED_TASKS}.')
        COMM.Abort()
    logging.debug('User configs checked; No aborts called')

    return


def prep_ligands():
    '''
    Function prep_ligands returns a list where each item is the path to a pickled
    and compressed text file containing multiple ligand strings (ignores pkl or
    dat format).

    Return: 
        ligand_paths (list): a list where each item containis multiple ligand strings
        (the items are pickled and compressed).
    '''
    
    ligand_paths = []
    for dirpath, _, filenames in os.walk(args.ligand_library):
        for filename in filenames:
            if filename.endswith('.pkl') or filename.endswith('.dat'):
                ligand_paths.append(f'{dirpath}/{filename}')
    logging.debug(f'Ligands prepped; Number of ligand batch files = {len(ligand_paths)}')

    return ligand_paths


def clean_as_we_go():
    '''
    This function is only called by Rank 1. It will continue to look at the results
    and remove bad-scoring molecules until it receives a 'stop working' signal from
    Rank 0.
    '''

    logging.debug('Entered the clean_as_we_go function')
    ligands_folder = 'output/results/ligands'
    file_counter = 0
    lines_to_keep = []

    while True:
        # Receive message from any source
        message = COMM.recv(source=MPI.ANY_SOURCE)

        # Check if it's a stop message from rank 0
        if message == 'stop working':
            logging.info('Rank 0 sent a message to rank 1 to stop')
            break

        # check if its a file result update
        elif 'File results' in message:
            file_counter += 1
            current_time = time.strftime('%H:%M:%S')
            sender_rank = message.split('_')[1]  
            logging.debug(f'{message} at {current_time}')
            if file_counter >= FILES_BEFORE_CLEANUP:
                logging.info(f'file counter hit {file_counter}; starting clean up')
                sort_for_rank1()
                lines_to_keep = delete_files(NUMBER_OF_OUTPUTS)
                file_counter = 0

    # Inform Rank 0 that Rank 1 is done 
    COMM.send('Stopping clean_as_we_go', dest=0, tag=1)
    COMM.send(lines_to_keep, dest=0, tag=2)
    logging.info('Rank 1 finished cleaning and sent completion message to Rank 0')

    return


def sort_for_rank1():
    '''
    This function cats all results files into one, it arranges ligands based on
    the highest score, the sorted results are written to sorted_scores_all.txt,
    finally cleans up the directory
    '''

    logging.info('entering sort_for_rank1')
    try:
        os.remove('./output/results/sorted_scores_all.txt')
    except:
        logging.warning('sorted_scores_all.txt does not exist')
    try:
        os.remove('merged_results.txt')
    except:
        logging.warning('merged_results.txt does not exist')

    subprocess.run(['cat results_*.txt >> merged_results.txt'], shell=True)
    inputfile = 'merged_results.txt'
    outputfile = './output/results/sorted_scores_all.txt'
    result = []

    with open(inputfile, 'r') as fin:
        line = fin.readline()    
        while line:
            filename = basename(line.split()[-1])
            try:
                v = fin.readline().split()[0]
            except IndexError as e:
                logging.debug('probably race condition')
                subprocess.run(['cat merged_results.txt'], shell=True)
            logging.debug(v)
            result.append(f'{v} {filename}\n')
            line = fin.readline()

    with open(outputfile, 'w') as fout:
        fout.writelines(sorted(result, key=lambda x: float(x.split()[1])))

    return


def delete_files(num):
    '''
    Open the file with all scores sorted and delete any ligands that are not 
    in the top N (assuming N = NUMBER_TO_KEEP and length of file > N)

    Args:
        num - number of docked ligands to keep

    Return:
        lines_to_keep - list of top N ligands
    '''

    with open('./output/results/sorted_scores_all.txt', 'r') as fin:
        lines = fin.readlines()
        logging.info(f'number of lines = {str(len(lines))}')
            
    # Remove Top N elements
    lines_to_remove = lines[num:]
    lines_to_keep = lines[:num]
    logging.info(f'number of lines in linestoremove={len(lines_to_remove)}')

    try:
       for dirpath, _, filenames in os.walk('./output/pdbqt'):
               for line in lines_to_remove:
                   name = line.split()[0]
                   for filename in filenames:
                       file_path = os.path.join(dirpath, filename)
                       if filename == f'output_{name}': # this is correct
                           os.remove(file_path)
                           logging.debug(f'file path {file_path} was removed')
    except:
        logging.warning('file could not be deleted')
                
    return lines_to_keep  


def processing():
    '''
    Function processing handles sending ligands to worker ranks, receiving ligands
    from rank 0, and running docking for each ligand set.
    '''

    if args.forcefield == 'vina':
        v = Vina(sf_name='vina', cpu=CPUS, verbosity=VERBOSITY)
        if FLEXIBLE == True:
            v.set_receptor(f'{RECEPTOR}.pdbqt', f'{FLEX_RECEPTOR}.pdbqt')
        else:
            v.set_receptor(f'{RECEPTOR}.pdbqt')
        v.compute_vina_maps(center=[CENTER_X, CENTER_Y, CENTER_Z],
                            box_size=[SIZE_X, SIZE_Y, SIZE_Z])
    elif args.forcefield == 'ad4':
        v = Vina(sf_name='ad4', cpu=CPUS, verbosity=VERBOSITY)
        v.load_maps(map_prefix_filename = RECEPTOR)
        
    # Ask rank 0 for ligands and dock until rank 0 says done
    count = 1
    directory = 1
    while True:
        COMM.send(RANK, dest=0) # Ask rank 0 for another set of ligands ###
        ligand_set_path = COMM.recv(source=0) # Wait for a response 

       # COMM.send(RANK,dest = 1)

        if ligand_set_path == 'no more ligands':
            current_time = time.strftime('%H:%M:%S')
            logging.debug(f'Rank {RANK} has finished all work at {current_time}')
            COMM.send('message received--proceed to post-processing',dest=0)
            break
        try:
            ligands = unpickle_and_decompress(ligand_set_path)
        except Exception as e:
            logging.error(f'Error on rank {RANK}: could not \
                            unpickle/decompress {ligand_set_path}.')
            logging.debug(e)
        try:
            run_docking(ligands, v, directory)
        except Exception as e:
            #logging.error(f'Error on rank {RANK}: docking error with \
            #                ligand set {ligand_set_path}, ligands {ligands}.')
            logging.error(f'Error on rank {RANK}: docking error with ligand set {ligand_set_path}.')
            logging.debug(e)

        count += 1
        if count == 100:
            count = 1
            directory += 1

    return


def unpickle_and_decompress(path_to_file):
    '''
    Given an input filepath, decompress and unpickle the file. Return the contents
    of the file as a dictionary where keys are ligand filenames and values are the
    actual ligand strings

    Args:
        path_to_file - file path

    Return:
        dictionary_of_ligands - dictionary of ligands
    '''

    with open(path_to_file, 'rb') as f:
        compressed_pickle = f.read()

    try:
        depressed_pickle = blosc.decompress(compressed_pickle)
    except Exception as e:
        logging.error(f'Error decompressing ligand batch {compressed_pickle}; Error = {e}')
        depressed_pickle = ''

    try:
        dictionary_of_ligands = pickle.loads(depressed_pickle)
    except Exception as e:
        logging.error(f'Error unpickling ligand batch {depressed_pickle}; Error = {e}')
        dictionary_of_ligands = {}

    return dictionary_of_ligands


def run_docking(ligands, v, directory):
    '''
    Run AutoDock on each liagand in the given set. Output a pdbqt file with the pose
    and all scores. Append the ligand name (filename) and its best pose/score to a 
    temporary results_N.txt file

    Args:
        ligands - dictionary of ligands
        v - Vina instance
        directory - path to file
    '''

    if not bool(ligands):
        logging.error('Ligand batch was set to empty; no ligands docked')
        return

    output_directory = f'./output/pdbqt/{RANK}{directory}'
    if not exists(output_directory):
        os.makedirs(output_directory)

    for _, filename in enumerate(ligands):
        ligand = ligands[filename]
        if not bool(ligand) or ligand == {}:
            logging.error(f'Ligand {filename} string was empty; skipping')
            continue
        try:
            v.set_ligand_from_string(ligand)
        except Exception as e:
            logging.error(f'Error setting ligand {filename}; Error = {e}')
            continue
        try:
            v.dock(exhaustiveness=EXHAUSTIVENESS)
        except Exception as e:
            logging.error(f'Error docking ligand {filename}; Error = {e}')
            continue
        try:
            v.write_poses(f'{output_directory}/output_{filename}',
                           n_poses=POSES, overwrite=True)
        except Exception as e:
            logging.error(f'Error writing ligand {filename}; Error = {e}')
            continue

        ### putting subprocess write function on one line to avoid race condition
        #subprocess.run([f"grep -i -m 1 'REMARK VINA RESULT:' \
        #               {output_directory}/output_{filename} \
        #               | awk '{{print $4}}' >> results_{RANK}.txt; echo {filename} \
        #               >> results_{RANK}.txt"], shell=True)
        command=f'grep -i -m 1 "REMARK VINA RESULT:" {output_directory}/output_{filename} \
                  | awk -v var="{filename}" "{{print \$4}} END{{print var}}" >> results_{RANK}.txt'
        subprocess.run([command], shell=True)

        COMM.send(f'File results_{RANK}.txt was generated', dest=1)

    return


def final_sort(start_time):
    '''
    The final_sort function cats all results files into one, it arranges ligands based
    on the highest score; prints these sorted results are written to sorted_scores.txt;
    finally cleans up the directory
    '''

    try:
        os.remove('./output/results/sorted_scores_all.txt')
    except:
        logging.warning('sorted_scores_all.txt does not exist')
    try:
        os.remove('./output/results/sorted_scores.txt')
    except:
        logging.warning('sorted_scores.txt does not exist')
    try:
        os.remove('merged_results.txt')
    except:
        logging.warning('merged_results.txt does not exist')
    try:
        os.remove('./output/results/total_time.txt')
    except:
        logging.warning('total_time.txt does not exist')

    subprocess.run(['cat results_*.txt >> merged_results.txt'], shell=True)
    inputfile = 'merged_results.txt'
    outputfile = './output/results/sorted_scores.txt'
    outputfile_all = './output/results/sorted_scores_all.txt'
    outputfile_time = './output/results/total_time.txt'
    result = []

    with open(inputfile, 'r') as fin:
        line = fin.readline()
        while line:
            filename = basename(line.split()[-1])
            v = fin.readline().split()[0]
            result.append(f'{v} {filename}\n')
            line = fin.readline()

    sorted_result = sorted(result[:], key=lambda x: float(x.split()[1]))

    with open(outputfile, 'w') as fout:
        fout.writelines(sorted_result[:NUMBER_OF_OUTPUTS])

    with open(outputfile_all, 'w') as fout:
        fout.writelines(sorted_result[:])

    end_comp_time = time.time()
    tot_time = end_comp_time - start_time
    with open(outputfile_time, 'w') as fout:
        fout.write(f"Total compute time: {tot_time}\n")

    return


def isolate_output():
    '''
    This function copies the user-specified top N ligand output files to a 
    single directory, and also cats all the top N ligand files into one big
    file, then tars the results
    '''

    top_ligand_filenames = []
    ligands_docked = 0
    logging.debug('running isolate_output')

    with open('./output/results/sorted_scores.txt', 'r') as fin:
        for _, line in enumerate(fin):
            top_ligand_filenames.append(line.split()[0])

    for dirpath, _, filenames in os.walk('./output/pdbqt'):
        ligands_docked += len(filenames)
        for top_filename in top_ligand_filenames:
            for filename in filenames:
                if filename == f'output_{top_filename}':
                    shutil.move(f'{dirpath}/{filename}', './output/results/ligands')
    
    fout = open('./output/results/combined_docked_ligands.pdbqt', 'w+')
    while top_ligand_filenames:
        with open(f'./output/results/ligands/output_{top_ligand_filenames.pop(0)}', 'r') as fin:
            lines = fin.read()
            fout.write(lines)
    fout.close()
    
    logging.debug(f'Total ligands saved = {ligands_docked}')
    subprocess.run([f'tar -czf results.tar.gz ./output/results'], shell=True)

    return


def reset():
    '''
    This function resets directory by removing all temp files from this job
    '''

    for dirpath, _, filenames in os.walk('.'):
        for filename in filenames:
            if filename.endswith(('.map', '.txt', '.gpf', '.fld', '.xyz')):
                logging.debug(f'Removed {filename}\n')
                os.remove(f'{dirpath}/{filename}') 
    shutil.rmtree('./output')
    shutil.rmtree('./configs')

    return


def main():
    '''
    Main loop:
        Rank 0: preprocess, manage work, postprocess
        Rank 1: clean as we go
        Rank 2-N: dock molecules
    '''

    current_responses = 0

    if RANK == 0:
        start_time = time.time()
        logging.info('Beginning main loop')

        # Pre-Processing
        pre_processing()
        ligands = prep_ligands()
        logging.info('Pre-processing finished; Starting SENDRECV with all ranks')
        
        # Let other ranks know pre-processing is finished; they can now ask for work
        for i in range(1, SIZE):
            COMM.send('pre-processing finished; ask for work', dest=i)
        logging.info('SEND finished; All ranks sent')
        
        pre_responses = 0
        while pre_responses != (SIZE-1):
            response = COMM.recv(source = MPI.ANY_SOURCE)
            if response == 'ready to go':
                pre_responses += 1
        logging.info('RECV finished; All ranks responded; Starting main work')

        COMM.send('Rank 1 will be ready to clean as we go', dest=1, tag=1)
        logging.info('Told rank 1 to clean as we go')

        # Until all ligands have been docked, send more work to worker ranks
        while ligands:
            source = COMM.recv(source=MPI.ANY_SOURCE)
            COMM.send(ligands.pop(), dest=source)
        logging.info('List of ligands is now empty')

        # When all ligands have been sent, let worker ranks know they can stop
        for i in range(2,SIZE):
                COMM.send('no more ligands', dest=i)
        logging.info('All ranks informed there is no more work')
        
        current_responses = 0
        while current_responses != (SIZE-2):
            response = COMM.recv(source = MPI.ANY_SOURCE)
            if response == 'message received--proceed to post-processing':
                current_responses += 1
        logging.info(f'Ranks 2-N have responded; Telling Rank 1 to stop working')
        
        COMM.send('stop working', dest=1)
        finished_message = COMM.recv(source=1, tag=1)
        logging.info(finished_message)
        top_ligands_message = COMM.recv(source=1, tag=2)
        logging.info(f'Rank 1 has responded; Proceeding to post-processing')

        # Post-Processing
        final_sort(start_time)
        isolate_output()
        reset()
        logging.info('Finished sorting and reset temp files')
        end_time = time.time()
        total_time = end_time - start_time

        logging.info(f'Script runtime = {total_time}.')
        logging.info(f'Nodes: {NODES}.')
        logging.info(f'Tasks: {TASKS}.')
        logging.info(f'Library: {LIBRARY_NAME}.')

    elif RANK == 1:
        COMM.recv(source=0) # Wait for rank 0 to finish pre-processing
        COMM.send('ready to go', dest=0)
        message = COMM.recv(source = 0)
       
        if message == 'Rank 1 will be ready to clean as we go':
            logging.info('Beginning clean as we go')
            clean_as_we_go()

    else: # All ranks besides rank 0
        COMM.recv(source=0) # Wait for rank 0 to finish pre-processing
        COMM.send('ready to go', dest=0)
        logging.info('Beginning processing')
        processing()    

    return


if __name__ == '__main__':
    main()


APPID ?= wallen-autodock-vina
APPVER ?= 1.2.3
# DOCKERID ?= wjallen/autodock-vina:1.2.3
# autodock-joey IS USED FOR TESTING. DO NOT DEPLOY
DOCKERID ?= joeleehen/autodock-joey:0.0.1

deploy:
	python3 ../utils/manage_apps app.json

build:
	docker build -t $(DOCKERID) .

run:
	docker run --rm -it --entrypoint=/bin/bash $(DOCKERID)

push:
	docker push $(DOCKERID)


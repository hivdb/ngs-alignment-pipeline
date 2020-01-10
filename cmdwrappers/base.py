#! /usr/bin/env python

import click
import docker

DOCKER_IMAGE = 'hivdb/five-prime-alignment-exp:latest'
REFINIT_FUNCTIONS = {}
ALIGN_FUNCTIONS = {}
AUTOREMOVE_CONTAINERS = False


def docker_execute(command, volumes):
    client = docker.from_env()
    try:
        return client.containers.run(
            DOCKER_IMAGE,
            command=command,
            volumes=volumes,
            remove=AUTOREMOVE_CONTAINERS,
            stderr=True,
            detach=False)
    except docker.errors.ContainerError as ex:
        click.echo(
            ex.stderr.decode('U8'),
            err=True)
        exit(1)


def refinit_func(name):

    def wrapper(func):
        REFINIT_FUNCTIONS[name] = func
        return func

    return wrapper


def align_func(name):

    def wrapper(func):
        ALIGN_FUNCTIONS[name] = func
        return func

    return wrapper


def get_programs():
    return sorted(ALIGN_FUNCTIONS.keys())


def get_refinit(name):
    return REFINIT_FUNCTIONS[name]


def get_align(name):
    return ALIGN_FUNCTIONS[name]

# Intall 'pytest' command to run tests available in anywhere in the current directory or in the subdirectories.
import setuptools
from setuptools import setup
from setuptools.command.test import test as testcommand

class PyTest(testcommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.tests")]

    def initialize_options(self):
        testcommand.initialize_options(self)
        self.pytest_args = []

    def run_tests(self):
        import pytest
        import sys
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)


setup(cmdclass={'tests': PyTest})
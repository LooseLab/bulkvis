from setuptools import setup

version = {}
with open("bulkvis/_version.py") as fh:
    exec(fh.read(), version)

setup(
    name="bulkvis",
    version=version["__version__"],
    author="Alexander Payne",
    entry_points={
        'console_scripts': [
            'bulkvis=bulkvis.bulkvis:main',
        ],
    },
    package_dir={'bulkvis': 'bulkvis'},
)

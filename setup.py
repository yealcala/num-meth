from setuptools import setup, find_packages

setup(
    name="num-meth",
    version="0.4",
    author="Yeray",
    description="Implementations of 'Métodos Algorítmicos en Matemáticas'",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    install_requires=[
        # "sagemath>=10.6"
    ],
    license="MIT"
)

# python -m pip install -e .

#  python -m build
#  twine upload --verbose --skip-existing --repository testpypi dist/*
#  twine upload dist/*
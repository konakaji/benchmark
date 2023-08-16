import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="benchmark",
    version="0.0.2",
    author="kouhei nakaji",
    author_email="nakajijiji@gmail.com",
    description="You can receive the message 'Hello!!!'",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/konakaji/benchmark",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "qwrapper @ git+ssh://git@github.com/konakaji/qwrapper.git",
        "tequila-basic==1.8.9",
        "openfermion>=1.5.1",
        "openfermionpyscf>=0.5",
        "pyscf==2.0.1"
    ],
    python_requires='>=3.6',
)

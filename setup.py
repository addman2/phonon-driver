import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pd",
    version="0.0.1",
    author="Otto Kohulak",
    author_email="pravod@gmail.com",
    description="phonon-driver",
    url="https://gitlab.mtf.stuba.sk/kohulak/phonon-driver",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

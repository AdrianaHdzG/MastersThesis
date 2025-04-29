from setuptools import setup, find_packages
import os

ASREwinlib = os.path.join("ASREpy","ASREcpp", "bin", "win32", "ASRElib.dll")
ASREmacOSlib = os.path.join("ASREpy","ASREcpp", "bin", "macOS", "libASRElib.dylib")

setup(
    name="ASREpy",
    version="0.1.0",
    author="Jinyan Zhao",
    author_email="jinyan_zhao@berkeley.edu",
    description="A Python package for Analysis of Structural Response to Excavation (ASRE)",
    # long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/your-repo",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        # List your package dependencies here, e.g.,
        # "numpy >= 1.18.0",
    ],
    # scripts=[ASREwinlib],
    include_package_data=True,
    # package_data= {
    #    'ASREpy': ['ASREpy\\ASREcpp\\bin\\win32\\ASRElib.dll',
    #               'ASREpy\\ASREcpp\\bin\\macOS\\libASRElib.dylib']
    # },
    package_data= {
       'ASREpy': [ASREwinlib,
                  ASREmacOSlib]
    },
    
)

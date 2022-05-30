from setuptools import setup
from setuptools import find_packages

setup(

    name = "score3pHLA",
    version = "0.2.1",
    description = "3pHLA-score: rank your pHLA structures",
    author = "Anja Conev",
    author_email = "ac121@rice.edu",
    url = "https://github.com/KavrakiLab/3pHLA-score",
    long_description=open("README.md").read(),
    long_description_content_type='text/markdown',
    packages = find_packages(),
    package_dir = {"score3pHLA":'./score3pHLA/'},
    install_requires=["biopython == 1.79",
                    "scikit-learn>=0.22.2",
                    "numpy >= 1.19"],
    package_data={"score3pHLA": ["*.pdb", 
                                "./3pHLAmodels/*.pkl"]},
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "score3pHLA = score3pHLA.__main__:main"
        ]
    },
)
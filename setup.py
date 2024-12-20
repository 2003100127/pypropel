from setuptools import setup, find_packages

setup(
    name="pypropel",
    version="0.1.3",
    keywords=["conda", "pypropel"],
    description="processing protein data",
    long_description="processing protein data",
    license="GPL v3.0",

    url="https://github.com/2003100127",
    author="Jianfeng Sun",
    author_email="jianfeng.sun@ndorms.ox.ac.uk",

    packages=find_packages(),
    include_package_data=True,
    platforms="any",
    python_requires='>3.9',
    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        'seaborn',
        'biopython',
        'scikit-learn',
        # 'pyfiglet', # ==0.8.post1
    ],
    entry_points={
        'console_scripts': [
            'pypropel=pypropel.main:main',
        ],
    }
)
import setuptools

with open("README", "r") as readme:
    long_description = readme.read()


setuptools.setup(
        name="surfgraph",
        version="0.2.0",
        scripts=['bin/analyze_chem_env.py', 'bin/generate_sites.py'],
        author="Tristan Maxson, Siddharth Deshpande",
        author_email="tgmaxson@gmail.com, sidd20111992@gmail.com",
        description="A package for generating chemical environment graphs of surfaces",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://gitlab.com/jgreeley-group/graph-theory-surfaces",
        python_requires=">=3",
        install_requires=[
            'numpy>=1.11.3',
            'networkx>=2.5',
            'matplotlib>=2.0.0',
            'ase>=19.1',
        ],
        packages=setuptools.find_packages(),
        classifiers=[
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Chemistry",
            "Topic :: Scientific/Engineering :: Physics",
            "Topic :: Scientific/Engineering",
            'License :: OSI Approved :: '
            'GNU Lesser General Public License v2 or later (LGPLv2+)',
            "Framework :: Matplotlib"],
        )


import setuptools

with open("README", "r") as readme:
    long_description = readme.read()


setuptools.setup(
        name="surfgraph",
        version="0.1",
        scripts=['bin/analyze_chem_env.py', 'bin/generate_sites.py'],
        author="Tristan Maxson, Siddharth Deshpande",
        author_email="tgmaxson@gmail.com, sidd20111992@gmail.com",
        description="A package for generating chemical environment graphs of surfaces",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://gitlab.com/jgreeley-group/graph-theory-surfaces",
        packages=setuptools.find_packages(),
        classifiers=[
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Chemistry",
            "Topic :: Scientific/Engineering :: Physics",
            "Topic :: Scientific/Engineering",
            "Framework :: Matplotlib"],
        )


INSTALLATION PROCEDURE
======================

this instruction will guide you through the installation of all the libraries
that are needed to get a working copy of the LifeV library

1) clone the LifeV github repo

```bash
git clone https://github.com/lifev/lifev.git
```

1b) OPTIONAL - if you plan to use your own configuration script for lifev
you chould also clone the cmake repo. If you are going to use these scripts
it will be cloned automatically.

```bash
cd lifev
git clone https://github.com/lifev/cmake.git
```

2) move to the tools/install_scripts directory

```bash
cd tools/install_scripts
```

3) copy config_example in config and edit it for your needs

```bash
cp config_example config
edit config
```
this file contains the basic configurations for the installation, such as
the installation directory.

4) copy libpath_example in libpath and edit it for your needs

```bash
cp libpath_example libpath
edit libpath
```

this file sets the configuration for the libraries that will be used or
installed.

5) run the build script

```bash
./build
```


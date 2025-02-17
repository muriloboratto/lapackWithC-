# Sample C++ project using LAPACK

Based on [LAPACK codes](https://www.netlib.org/lapack/), this C++ project aims to:

* Compile on cross-platform software using make
* Publish software documentation to GitHub Pages https://muriloboratto.github.io/cpp-example/


## Building

> :information_source: These instructions describe how to run CMake commands manually. 
1. Install `git`, a C++ compiler, `CMake`, and a build system available as a `CMake` generator
2. Download this repository.
    ```bash
    $ git clone https://github.com/muriloboratto/lapackWithCplusplus
    $ cd lapackWithC++
    ```
3. Run the configuration step in the build directory. 
    ```bash
    $ make
    ```
    
## Running

After being built, the binaries can be found in the build directory:
```bash
$ ./principal <size problem>
$ ./principal     10
``

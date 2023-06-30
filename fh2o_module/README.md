# Build

This module must be build in every system.

Install dependences

```bash
apt install python3-numpy gfortran
```

Run the following command

```bash
python3 -m numpy.f2py -c  li_dawes_guo.f fh2o-pipnn.f -m li_dawes_guo
```

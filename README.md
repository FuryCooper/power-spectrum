### power spectrum

```
git clone https://github.com/FuryCooper/power-spectrum.git
```

### Installation
Before doing "make", you might need to modify Makefile to include your own fftw library. Also kindly remind you
that fftw3 might be not supported, so fftw2 library is recommended.

```
make
make install
```

### how to run it 
```
./POWERSPEC 4 ./example/example.param
```

4 is your thread number. you can set it individually.

### Uninstallation

```
make clean
```
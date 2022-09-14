# Root Function

The idea is simple, receive a generic function and return the position of the
root using the newtons method.

> What's the newton method?

Check the following gif to understand in an intuitive way.

![newton-method](https://upload.wikimedia.org/wikipedia/commons/thumb/e/e0/NewtonIteration_Ani.gif/300px-NewtonIteration_Ani.gif)

To undestand in detail read the wikipedia page [here](https://en.wikipedia.org/wiki/Newton%27s_method)

## Benchmark

Just a simples bech with newtons_method 2 write in `python` and `golang`

```bash
# Python
time python3 ./newtons-method-2.py 
Converge in (0.6234994049308765, 0.02803775852868567)

real 0m0,024s
user 0m0,024s
sys  0m0,000s
```

```bash
# Golang
time ./newtons-method-2 
Converge in [0.6234994049308766 0.028037758528685664] 

real 0m0,003s
user 0m0,000s
sys  0m0,003s
```

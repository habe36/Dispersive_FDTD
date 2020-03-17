## Lorentz_Dispersive_FDTD
Finite Difference Time Domain Method for solving Lorentz Dispersive Media with TGM and ADE Method.
These are Python scripts (*.py) and GNUPlot Scripts (*.gp).
lorentz1d.py and tgm1d.py are solvers to solves one dimensional electromagnetic wave reflection with Lorentz Dispersive Media.
lorentz1d is coded with ADE method, a convensional method and tgm1d is coded with Transient Green Method.
## Usage
> `% python tgm1d.py 20`

to solve numerically for 20GHz pole Lorentz media. After a while this generate a file, "tgm1d.data".
And then,
>% python lorentz1d.py 20

to generate "lorentz1d.data"
In order to calculate the reflection coefficients, use reflect.py.
> % python reflect.py tgm1d.data

> % python reflect.py lorentz1d.data

Use gnuplot to plot the result.
> % gnuplot

> load 'reflection.gp'

to display the result.
## License
GPL3.

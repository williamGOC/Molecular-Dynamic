# Molecular Dynamic
This code can be used for the study of classical gases with weak interaction (Van der Waals). The evolution of the system is done by the ```Speed Verlet algorithm```.
#### Compilation

The compilation of the codes ```c ``` is very simple, we simply use:

```
gcc input.c -o output.out -lm
./output.out
```
#### Pipe in gnuplot

In ```linux``` you can call ```gnuplot``` from your ```c``` code as follows:

```c
FILE *pipe = popen("gnuplot", "w");
assert(pipe);
/*
Development code of our program
*/
pclose(pipe);
```
This form of call can be used in a ```window```, provided that the compilation of the code is done through a console (```gnuplot``` and ```gcc``` must be environment variables).

```c
/*
Print that allows to visualize the evolution of the system.
Keyword Arguments:
sys: class_system * (particle system).
gnuplot: ----> environment variable (Window).
         ----> object FILE * (Linux).
*/
void printer_system(class_system *sys, FILE *gnuplot){
    int i;
   
    fprintf(gnuplot, "set title '{/=20 Gas Model, pass %d}'\n", sys->accountant);
    fprintf(gnuplot, "set xlabel  '{/=15 X}'\n");  
    fprintf(gnuplot, "set zlabel  '{/=15 Z}'\n");
    fprintf(gnuplot, "set ylabel  '{/=15 Y}'\n");
    fprintf(gnuplot, "splot '-' with p pointtype 6 t ''\n");
  
    for(i = 0; i < N; i++)
        fprintf(gnuplot, "%lf\t%lf\t%lf\n", sys -> position[i][0], 
                                    sys -> position[i][1], sys ->position[i][2]);
    
    fprintf(gnuplot, "e\n");
    fflush(gnuplot);
}
```

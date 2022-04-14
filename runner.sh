gcc -o poisson poisson.c -pthread
for ii in {1..1}
do
   ./poisson -n 301 -i 500 >out.txt
done

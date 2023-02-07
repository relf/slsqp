To generate compile commands json file required by c2rust in top directory

> mkdir build
> cd build
> cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 ../slsqp
> cp compile_commands.json ../slsqp/compile_commands.json

Then rust generation from C

> cd ../slsqp
> c2rust transpile compile_commands.json 
# treeCmp

TreeCmp: comparison of trees in polynomial time. See the [manual](TreeCmp_manual.pdf) for usage.

## Compilation from the command line

Normally, compilation is done via an IDE like NetBeans or Eclipse

To compile via the command-line, creating files in the directory `out/class`, you should be able to do the 
following from the top level folder

```
javac -d out/class -cp lib/commons-cli-1.4.jar src/treecmp/*.java src/treecmp/*/*.java src/pal/*/*.java
```

The resulting compiled files can be run directly, for example by issuing the command
`java -cp out/class:lib/commons-cli-1.4.jar treecmp.Main`
(you will need to replace the colon with a semicolon on windows systems). However, it is usually easier to create a stand-alone .jar executable as below.

## Creating a jar executable from the command line

The jar executable (e.g. `TreeCmp.jar`) should be created in the `bin` directory, as it expects to find the `config/config.xml` file one level above the directory in which the executable resides. Once the class files have been compiled in the `out/class` directory as above, this jar file can be created in the correct place using the MANIFEST.MF file in `src/META-INF`:

```
jar cvfm bin/TreeCmp.jar src/META-INF/MANIFEST.MF -C out/class/ .
```

As defined in the manifest, extra libraries are expected to be placed in `lib` directory in the same place as the jar file. An easy way to do this is to move or copy the `lib` folder (containing `commons-cli-1.4.jar`) into the `bin` directory. For example, on unix-like systems you could do `cp -a lib bin/`, on windows `xcopy lib bin\`.

Once the `commons-cli-1.4.jar` library has been copied to the correct place, the `TreeCmp.jar` file can be run as described in the [manual](TreeCmp_manual.pdf) (e.g. 
`java -jar bin/TreeCmp.jar -w 2 -d ms -i examples/beast/testBSP.newick -o testBSP.newick_w_2.out -I`

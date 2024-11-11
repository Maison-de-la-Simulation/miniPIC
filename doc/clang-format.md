# Clang-format

`clang-format` is a command line tool, part of the clang set of tools, that allows
developers to format C/C++ files into common or personalized styles.

This code is formatted following the rules provided by the file `.clang-format`.

You can easily lear how to use `clang-format` using the foollowing links:

- [How to use clang-format alone or as a plugin with your favorite tools](https://clang.llvm.org/docs/ClangFormat.html)
- [How to tune the YAML format file](https://clang.llvm.org/docs/ClangFormatStyleOptions.html).

How to use clang-format with git is explained [on this link](https://electronjs.org/docs/development/clang-format)

### Formatting a file manually

The file `.clang-format` has to be in the parent file of the file
or where you run the command.

To test the formatting process on the `.cpp` or `.h` file of your choice:
```bash
clang-format -style=file /src/<file of your choice>
```

The result is shown in your terminal.

In order to directly change the file without terminal output, use the flag `-i`:
```
clang-format -style=file -i /src/<file of your choice>
```

### Compatibility with git

In order to format automatically the modified files seen by git (before doing a commit for instance), you can run the command
```
git-clang-format
```

### Reformatting all the sources

Cpp files:

```bash
find src/* -name \*.cpp -exec clang-format -style=file -i '{}' \;
```

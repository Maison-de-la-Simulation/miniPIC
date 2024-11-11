# CI

The code is tested using the validation processes described in [the Validation section](./validation.md)

Code validation is automatically performed on `master` and `develop` branches and pull requests.

For other branches, you can trigger the CI manually putting `[run ci]` in the commit message.

You can as well disable the CI by putting `[skip ci]` or [`draft`] in the commit message.
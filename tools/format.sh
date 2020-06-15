#! /usr/bin/env bash
find . -regex '.*\.\(cpp\|hpp\|cc\|cxx\|h\)' -exec clang-format-9 -style=file -i {} \;
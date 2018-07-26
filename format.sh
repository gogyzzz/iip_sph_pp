find ./ -iname *.c -o -iname *.h -o -iname *.tpp -o -iname *.cpp | xargs clang-format -i

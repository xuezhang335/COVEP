from deephlapan_main import *
from parse_args import CommandLineParser


def main():
    (opt, _) = CommandLineParser()
    deephlapan_main(opt)


if __name__ == '__main__':
    main()

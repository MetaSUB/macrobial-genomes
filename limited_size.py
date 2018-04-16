from sys import stdin
import click


def parseFasta(baseLim):
    curLen = 0
    for line in stdin:
        if curLen > baseLim:
            break
        line = line.strip()
        if len(line) == 0:
            continue
        elif line[0] != '>':
            curLen += len(line)
        yield line


@click.command()
@click.argument('base_lim', type=int)
def main(base_lim):
    for line in parseFasta(base_lim):
        print(line)

if __name__ == '__main__':
    main()

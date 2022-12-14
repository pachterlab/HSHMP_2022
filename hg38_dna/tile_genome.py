#!/usr/bin/env python3

import argparse
import gzip

def write_entry(out, head, body, FA_WIDTH=80):
    out.write(str.encode(f'{head}\n'))
    # Break down the body into chunks of length FA_WIDTH for legibility
    out.write(str.encode('\n'.join(body[i:i+FA_WIDTH]
                         for i in range(0, len(body), FA_WIDTH))))
    out.write(str.encode('\n'))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create a tiling cover of the\
                                                  supplied assembly scaffold file.')
    parser.add_argument('--length', type=int, default=500, help='Length of each tile')
    parser.add_argument('--fa', type=str, help='Genome assembly')
    parser.add_argument('--out', type=str, help='Output file')

    args = parser.parse_args()

    out = gzip.open(args.out, 'wb')

    with gzip.open(args.fa, 'rb') as fh:

        tile = ''
        tile_n = 1
        gene = ''
        for line in fh:

            l = line.decode('utf-8').rstrip('\n')
            if l[0] == '>':
                # Start new scaffold

                if gene != '' and len(tile) > 0:
                    # Write the last tile in last scaffold if this is not the
                    # first scaffold and there is a tile to be written
                    write_entry(out, f'{gene}.{tile_n}', tile)

                gene = l.split()[0]
                tile_n = 1
                tile = ''
                continue

            i = 0
            # Spool past Ns
            while i < (len(l)-1) and l[i] == 'N':
                i += 1
            if i == len(l) - 1:
                continue
            l = l[i:]

            left = args.length - len(tile)
            endp = min(left, len(l))
            tile += l[:endp]
            if len(tile) == args.length:
                write_entry(out, f'{gene}.{tile_n}', tile)
                tile_n += 1
                tile = l[endp:]

        # Write last entry if length of tile was not an integer multiple of the
        # total length of the last scaffold
        if len(tile) > 0:
            write_entry(out, f'{gene}.{tile_n}', tile)
    out.close()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-n", type=int, help="Provide a number please",
                    nargs='?', default=0, const=0, dest='n')
parser.add_argument("-vidName", type=str, help="Provide a number please",
                    nargs='?', default='vid_')
parser.add_argument('--foo', action='store_const', const=2, default=42)
parser.add_argument("-framewindow", type=int, help="Provide the frame window, please",
                    nargs='+', default=[0,1])
args = parser.parse_args()
print(args.vidName)
print(args)

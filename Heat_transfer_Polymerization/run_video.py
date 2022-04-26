import glob, os, sys
import re
import argparse, subprocess
import moviepy.video.io.ImageSequenceClip
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def sort_names(image_files):
    file_names = [os.path.basename(string) for string in image_files]
    prefix = os.path.dirname(image_files[0])
    times = [re.findall("\d+\.\d+", string)[0] for string in file_names]
    pos = image_files[0].find(times[0])
    prefix = image_files[0][:pos]
    postfix = image_files[0][(pos + len(times[0])):]
    times = sorted([float (t) for t in times])
    # print(" ".join(times))
    print(f"Time: {times[0]}  {times[-1]}")
    image_files = [prefix + f"{t:.6f}" + postfix for t in times]
    return image_files

def make_video(image_folder, videoTime):
    basename = os.path.basename(image_folder)
    image_files = [os.path.join(image_folder,img)
                   for img in os.listdir(image_folder)
                   if img.endswith(".png")]
    fps = int(len(image_files)/videoTime)
    print("FPS:", fps)
    image_files = sort_names(image_files)
    try:
        clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
        clip.write_videofile(f'{image_folder}/../{basename}.mp4')
    except:
        print(f"ERROR in video generation: image_files: {image_files}")


def make_zip(tree_path, zip_path, mode='w', skip_empty_dir=False):
    print(f"Make folder: {zip_path}")
    with ZipFile(zip_path, mode=mode, compression=ZIP_DEFLATED) as zf:
        paths = [Path(tree_path)]
        while paths:
            p = paths.pop()
            if p.is_dir():
                paths.extend(p.iterdir())
                if skip_empty_dir:
                    continue
            zf.write(p)
parser = argparse.ArgumentParser()
parser.add_argument("-dirPattern", type=str, help="Provide the pattern of search directory, please",
                    nargs='?', default='*/')
parser.add_argument("-genPic", type=str2bool, nargs='?',
                    const=True, default=True,
                    help="Allow generation of pictures, please")
# parser.add_argument("-genPic", type=bool, help="Allow generation of pictures, please",
#                     nargs='?', default=True)
# parser.add_argument("-genVid", type=bool, help="Allow generation of video, please",
#                     default=True)
parser.add_argument("-genVid", type=str2bool, nargs='?',
                    const=True, default=True,
                    help="Allow remove of directories, please")
# parser.add_argument("-genZip", type=bool, help="Allow generation of Zip, please",
#                     default=True)
parser.add_argument("-genZip", type=str2bool, nargs='?',
                    const=True, default=True,
                    help="Allow remove of directories, please")
# parser.add_argument("-normDir", type=bool, help="Allow remove of directories, please",
#                     default=True)
parser.add_argument("-rmDir", type=str2bool, nargs='?',
                    const=True, default=True,
                    help="Allow remove of directories, please")
parser.add_argument("-time", type=float, help="Provide the frames per second, please",
                    nargs='?', default=10)
parser.add_argument("-fn", type=str, help="Provide the pvd file names without extension, please",
                    nargs='?', default='output_res')

args = parser.parse_args()
print(args)
dirPattern = args.dirPattern
videoTime = args.time
output_res = args.fn
genPic = args.genPic
genVid = args.genVid
genZip = args.genZip
rmDir = args.rmDir

maxTime = 10000000
dirnames = glob.glob(f'./{dirPattern}/', recursive = False)
dirnames2 = []
for dirn in dirnames:
    if glob.glob(f'./{dirn}/{output_res}.pvd', recursive = False):
        dirnames2.append(dirn)
print('Found pvd files in:',dirnames2)


for dirn in dirnames2:
    print('Dirname:', dirn)
    if genPic:
        command = f"cd {dirn}/ && (ln -s ../one_by_one.py . || echo 'Note: one_by_one.py link exists!') && pvpython one_by_one.py -filenames {output_res}.pvd >out 2>err"
        print(command)
        subprocess.run(command, capture_output=True, shell=True)
    for fieldName in ['u', 'T', 'alpha_doc', 'mu_cell']:
        dirname = f'{dirn}/pone/{fieldName}'
        if genVid:
            make_video(dirname, videoTime)
        if genZip:
            make_zip(dirname, f"{dirname}/../{fieldName}.zip", mode='w', skip_empty_dir=False)
        if rmDir:
            subprocess.run(f"rm -r {dirname}/", capture_output=True, shell=True)
# subprocess.run(f"echo "FINISHED..."; sleep {maxTime};", capture_output=True, shell=True)

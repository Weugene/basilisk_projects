import glob
filenames='*ux.png'
for file in glob.glob(filenames):
    print(file)

import ffmpeg
(
    ffmpeg
    .input(filenames, pattern_type='glob', framerate=10)
    .output('Mux.avi')
    .run()
)

filenames='*noLambda2.png'
(
    ffmpeg
    .input(filenames, pattern_type='glob', framerate=10)
    .output('MnoLambda2.avi')
    .run()
)

filenames='*tracer.png'
(
    ffmpeg
    .input(filenames, pattern_type='glob', framerate=10)
    .output('Mtracer.avi')
    .run()
)

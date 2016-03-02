#!/bin/bash

# more sophisticated script for rendering a movie that takes in dynamic
# camera positions and targets 

usage()
{
  echo "render_movie.sh <frame rate> <path to Moby render> <camera position file> <path to osg files/output> <camera target file> <movie file>"
  exit 
}

main()
{


  # commands like this can be used to modify the osg files
  #rpl -q "ColorMode DIFFUSE" "ColorMode AMBIENT_AND_DIFFUSE" $1/*.osg

  # renders the osg files to images in order expected by ffmpeg
  
  a=1
  for i in $DATA_PATH/driver.out-*.osg; do
    new=$(printf "img%04d.png" ${a})
    pnt=$(awk 'NR == n' n=${a} $CAMERA_FILE)
    tgt=$(awk 'NR == n' n=${a} $TARGET_FILE)
    [ -f $DATA_PATH/$new ] && echo "Already created $new" || screen -d -m $RENDER_PATH/moby-render $i -p $pnt -t $tgt $DATA_PATH/${new}
    let a=a+1
  done

  read -rsp $'Press any key to continue to making video... \n' -n1 key

  # render at 25fps 
  ffmpeg -r 100 -i $DATA_PATH/img%04d.png -f mp4 -q:v 0 -vcodec mpeg4 $MOVIE_FILE 
}

# check for proper number of arguments
[ "$#" -ne 5 ] && ( usage && exit)

echo "Input count: " $# " , Values: " $@
echo "<path to Moby render>" $1
echo "<path to osg files/output>" $2
echo "<camera position file>" $3
echo "<camera target file>" $4
echo "<movie file>" $5

RENDER_PATH=$1
DATA_PATH=$2
CAMERA_FILE=$3
TARGET_FILE=$4
MOVIE_FILE=$5

main

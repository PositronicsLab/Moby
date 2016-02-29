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
  for i in $4/driver.out-*.osg; do
    new=$(printf "img%04d.png" ${a});
    let a=a+1;
    #render-osg $i -p 0.55 -0.75 0.5 -t $(awk 'NR == n' n=${a} $1/com.mat) -s=$1/scene.osg $1/${new};
    $2/moby-render $i -p $3 $4 $5 -t $(awk 'NR == n' n=$(($6 * ${a})) $3) $6/${new};
  done

  # render at desired frame rate 
  ffmpeg -r $1 -i $6/img%04d.png -f mp4 -q:v 0 -vcodec mpeg4 $6
}

# check for proper number of arguments
[ "$#" -ne 6 ] && ( usage && exit) || main $1 $2 $3 $4 $5 $6 


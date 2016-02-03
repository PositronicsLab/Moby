#!/bin/bash

# simple script for rendering a movie that takes in static camera positions
# and targets

usage()
{
  echo "render_movie_simple.sh <frame rate> <path to Moby render> <camera position> <path to osg files/output> <camera target> <movie file>"
  echo " -- sample camera position: 0 5 10"
  echo " -- sample camera target: 0 0 0"
  exit 
}

main ()
{
  # ffmpeg requires particular numbering
  a=1
  for i in $6/driver.out-*.osg; do
  new=$(printf "img%04d.png" ${a});
  let a=a+1;
  #render-osg $i -p 0.55 -0.75 0.5 -t $(awk 'NR == n' n=${a} $1/com.mat) -s=$1/scene.osg $1/${new};
echo $2/moby-render $i -p $3 $4 $5 -t $7 $8 $9 $6/${new} 
  $2/moby-render $i -p $3 $4 $5 -t $7 $8 $9 $6/${new};
done

# number all images img0001.png, etc...
echo ffmpeg -r $1 -i $6/img%04d.png -f mp4 -q:v 0 -vcodec mpeg4 ${10}
ffmpeg -r $1 -i $6/img%04d.png -f mp4 -q:v 0 -vcodec mpeg4 ${10}
}


# check for proper number of arguments
[ "$#" -ne 10 ] && ( usage && exit) || main $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}



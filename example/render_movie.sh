find $1 -maxdepth 1 -name "img*.png" -type f -print | xargs rm

# this code parses Links data to produce a target location
grep "CoM_x = " $1/out.log > $1/com.mat
rpl -q "]';" "" $1/com.mat
rpl -q "CoM_x = [" "" $1/com.mat

#rpl -q "ColorMode DIFFUSE" "ColorMode AMBIENT_AND_DIFFUSE" $1/*.osg

a=1
for i in $1/driver.out-*.osg; do
  new=$(printf "img%04d.png" ${a});
  let a=a+1;
  #render-osg $i -p 0.55 -0.75 0.5 -t $(awk 'NR == n' n=${a} $1/com.mat) -s=$1/scene.osg $1/${new};
  render-osg $i -p 2 3 1  -t $(awk 'NR == n' n=$(($2 * ${a})) $1/com.mat) $1/${new};
done

# number all images img0001.png, etc...
ffmpeg -r 25 -i $1/img%04d.png -f mp4 -q:v 0 -vcodec mpeg4 $1/movie.mp4

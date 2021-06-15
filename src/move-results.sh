dest=$1
mkdir "$dest"

mv output "$dest/"
mv post "$dest/"

mkdir output
mkdir output/velocity
mkdir post

# make sure github doesnt remove these files
touch output/velocity/.gitignore
touch post/.gitignore

mv *.gp "$dest/"

cp sample_inp.in "$dest"

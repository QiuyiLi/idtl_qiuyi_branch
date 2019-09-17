counter=1
while [ $counter -le 30 ]
do
python3 main.py -i data/tree_sample2.txt
((counter++))
done

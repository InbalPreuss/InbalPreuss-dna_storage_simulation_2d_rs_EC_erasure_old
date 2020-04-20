# Run on AWS
1. Create AWS account.
2. Create an EC2 instance (follow the instructions)
3. Click in AWS on **Connect to instance** and follow the instructions (getting key)
4. Open terminal/Powershell
5. Open Powershell
```console
ssh -i "inbal.pem" ubuntu@ec2-18-221-210-214.us-east-2.compute.amazonaws.com
```

5. install pip:
```console
sudo apt update
sudo apt install python3-pip
```

7. install the project
```console
git clone https://github.com/inbalpreuss/DnaStorage
cd DnaStorage/
pip3 install -r requirements.txt
```

8. create input file
```console
mkdir -p data/testing
touch data/testing/input_text.dna
echo "inbal preuss" > data/testing/input_text.dna
```

9. run dna storage
```console
python3 -m dna_storage.main
cat data/testing/simulation_data.9.text_results_file.dna
```

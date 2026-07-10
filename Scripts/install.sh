#!/bin/bash
echo "Installing  dependencies..."
sudo apt update && sudo apt install -y python3-venv && sudo apt install scons -y && sudo apt install libboost-dev -y && sudo apt install libhdf5-de -y && sudo apt install python3-venv -y && sudo apt install doxygen -y

echo "Creating Project directory..."
mkdir cantera_wall
cd cantera_wall
echo "Creating virtual environment 'wall_env'..."
python3 -m venv can_env
echo "Cloning repository into 'dummy'..."
git clone --recursive https://github.com/JacobDannecker/cantera_ma.git cantera_ma
echo "Activating virtual environment..."
source can_env/bin/activate
echo "Installing requirements..."
pip3 install -r  cantera_ma/Scripts/requirements.txt
echo "Adding alias to ~/.bashrc..."
PROJECT_ABS_PATH="$(pwd)"
echo "alias wall_env='cd $PROJECT_ABS_PATH && source can_env/bin/activate'" >> ~/.bashrc
source ~/.bashrc
echo "Running scons build..."
cd cantera_ma 
echo "can_env/*" >> .gitignore
scons build prefix=$VIRTUAL_ENV -j$(nproc)
scons install prefix=$VIRTUAL_ENV -j$(nproc)
echo "Setup complete. Use can_env to change into your working directory and activate the virtual python environment."
echo "Use deactivate to deactivate the virtual environment."



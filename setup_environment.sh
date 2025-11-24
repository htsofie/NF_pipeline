#!/bin/bash

# Environment setup script for Nextflow Phosphorylation Analysis Pipeline
# This script sets up a Python virtual environment with all required packages

set -e  # Exit on any error

echo "Nextflow Phosphorylation Pipeline Setup"


# Check if Python 3 is available
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in PATH"
    echo "Please install Python 3.8 or higher"
    exit 1
fi

# Check Python version
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}' | cut -d'.' -f1,2)
PYTHON_MAJOR=$(echo $PYTHON_VERSION | cut -d'.' -f1)
PYTHON_MINOR=$(echo $PYTHON_VERSION | cut -d'.' -f2)

if [ "$PYTHON_MAJOR" -lt 3 ] || ([ "$PYTHON_MAJOR" -eq 3 ] && [ "$PYTHON_MINOR" -lt 8 ]); then
    echo "Error: Python 3.8 or higher is required. Found: $PYTHON_VERSION"
    exit 1
fi

echo "Python version: $(python3 --version)"

# Detect operating system
if [[ "$OSTYPE" == "darwin"* ]]; then
    OS="macos"
    echo "Detected OS: macOS"
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    OS="linux"
    echo "Detected OS: Linux"
else
    echo "Warning: Unknown OS type: $OSTYPE"
    OS="unknown"
fi

# Create virtual environment
echo ""
echo "Creating virtual environment..."
# Check if venv module is available first
if ! python3 -c "import venv" 2>/dev/null; then
        echo "Error: Python venv module is not available."
        if [[ "$OS" == "linux" ]]; then
            PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}' | cut -d'.' -f1,2)
            echo ""
            echo "On Debian/Ubuntu systems, you need to install the python3-venv package:"
            echo "  sudo apt install python${PYTHON_VERSION}-venv"
            echo ""
            echo "After installing, run this script again."
        elif [[ "$OS" == "macos" ]]; then
            echo ""
            echo "On macOS, ensure Python development tools are installed."
            echo "You may need to install Xcode command line tools:"
            echo "  xcode-select --install"
        fi
    exit 1
fi

# Check if venv exists and is valid
if [ -d "venv" ]; then
    if [ -f "venv/bin/activate" ] || [ -f "venv/Scripts/activate" ]; then
        echo "Virtual environment already exists. Skipping creation."
    else
        echo "Broken virtual environment detected. Removing and recreating..."
        rm -rf venv
        if ! python3 -m venv venv 2>&1; then
            echo ""
            echo "Error: Failed to recreate virtual environment."
            if [[ "$OS" == "linux" ]]; then
                PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}' | cut -d'.' -f1,2)
                echo "On Debian/Ubuntu systems, you may need to install python3-venv:"
                echo "  sudo apt install python${PYTHON_VERSION}-venv"
            fi
            exit 1
        fi
        echo "Virtual environment recreated successfully."
    fi
else
    # Create new venv (venv module check already done above)
    if ! python3 -m venv venv 2>&1; then
        echo ""
        echo "Error: Failed to create virtual environment."
        if [[ "$OS" == "linux" ]]; then
            PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}' | cut -d'.' -f1,2)
            echo "On Debian/Ubuntu systems, you may need to install python3-venv:"
            echo "  sudo apt install python${PYTHON_VERSION}-venv"
        elif [[ "$OS" == "macos" ]]; then
            echo "On macOS, ensure Python development tools are installed."
        fi
        exit 1
    fi
    echo "Virtual environment created successfully."
fi

# Activate virtual environment
echo ""
echo "Activating virtual environment..."
if [[ "$OSTYPE" == "darwin"* ]] || [[ "$OSTYPE" == "linux-gnu"* ]]; then
    source venv/bin/activate
else
    source venv/Scripts/activate
fi
PYTHON_CMD="python"
PIP_CMD="pip"

# Upgrade pip
echo ""
echo "Upgrading pip..."
$PIP_CMD install --upgrade pip setuptools wheel

# Install requirements
echo ""
echo "Installing Python packages from requirements.txt..."
$PIP_CMD install -r requirements.txt

# Test installation
echo ""
echo "Testing installation..."
$PYTHON_CMD -c "
import sys
try:
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import yaml
    import Bio
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner
    from Bio.Blast import NCBIXML
    print('✓ All packages imported successfully!')
    print(f'✓ Python version: {sys.version}')
    print(f'✓ pandas version: {pd.__version__}')
    print(f'✓ numpy version: {np.__version__}')
    print(f'✓ matplotlib version: {plt.matplotlib.__version__}')
    print(f'✓ biopython version: {Bio.__version__}')
    print('✓ Environment setup complete!')
except ImportError as e:
    print(f'✗ Import error: {e}')
    sys.exit(1)
"

# Check for BLAST+
echo ""
echo "Checking for BLAST+ installation..."
if command -v makeblastdb &> /dev/null && command -v blastp &> /dev/null; then
    echo "✓ BLAST+ is installed"
    makeblastdb -version | head -n 1
else
    echo "⚠ Warning: BLAST+ not found in PATH"
    echo "  Please install BLAST+ to use this pipeline:"
    echo "  - macOS: brew install blast"
    echo "  - Ubuntu/Debian: sudo apt-get install ncbi-blast+"
    echo "  - Or download from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download"
fi

echo ""
echo "Setup complete!"
echo ""
echo "To use this environment:"
if [[ "$OSTYPE" == "darwin"* ]] || [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo "  source venv/bin/activate"
else
    echo "  source venv/Scripts/activate"
fi
echo ""
echo "To run the pipeline:"
echo "  ./nextflow run main.nf"
echo ""
echo "To deactivate the environment:"
echo "  deactivate"
echo ""


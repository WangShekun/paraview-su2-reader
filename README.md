# SU2 Mesh Reader (ParaView Python Plugin) 🧩

A lightweight ParaView Python reader plugin for visualizing **SU2 mesh files (`.su2`)**.  
Shared for convenience: **no support, no warranty** ⚠️

This reader provides two output ports:
- **Port 0: Volume Mesh** (`vtkUnstructuredGrid`) 
- **Port 1: Boundary Meshes (grouped by marker)** (`vtkMultiBlockDataSet`, one block per SU2 marker) 

**Tested on ✅**
- Ubuntu 24.04 LTS + ParaView 6.0.1  
- Windows 11 + ParaView 5.13.3  

---

## Installation 🔧

1. Open ParaView  
2. `Tools → Manage Plugins...`  
3. `Load New...` and select `SU2Reader.py`  
4. (Recommended) enable **Auto Load** so it loads automatically on startup  

---

## Usage ▶️

1. `File → Open` and choose a `*.su2` file  
2. Click `Apply`  

After opening, you will have two outputs to visualize 👀:
- **Port 0: Volume mesh** (the whole mesh) 
- **Port 1: Boundary meshes grouped by marker** (MultiBlock; view/extract markers individually) 

---

## Disclaimer (Important) ⚠️

This project is shared for reference and convenience only and is provided **"AS IS"**.  
No support or warranty is provided. The author is not responsible for any issues, damages, or disputes arising from the use of this code.  
Using this project means you understand and accept these terms.

---

## License 📄

MIT License (see `LICENSE` in the repository root).

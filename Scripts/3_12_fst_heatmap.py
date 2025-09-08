#!/bin/bash
# function: plot the fst matrix
import numpy as np
import matplotlib.pyplot as plt

"""
the matrix
	N1	N2	N3	N4	S
N1	0	0.0220845818406098	0.0275283510049481	0.0180433223898367	0.0906391542748899
N2	0.0220845818406098	0	0.0214320666421353	0.0117360753509322	0.0707225123213458
N3	0.0275283510049481	0.0214320666421353	0	0.0178295099593012	0.0894257678938331
N4	0.0180433223898367	0.0117360753509322	0.0178295099593012	0	0.0691401820113223
S	0.0906391542748899	0.0707225123213458	0.0894257678938331	0.0691401820113223	0
"""

fig, ax = plt.subplots(figsize=(5,6))
matrix = np.array([[0, 0.0220845818406098, 0.0275283510049481, 0.0180433223898367, 0.0906391542748899],
                   [0.0220845818406098, 0, 0.0214320666421353, 0.0117360753509322, 0.0707225123213458],
                   [0.0275283510049481, 0.0214320666421353, 0, 0.0178295099593012, 0.0894257678938331],
                   [0.0180433223898367, 0.0117360753509322, 0.0178295099593012, 0, 0.0691401820113223],
                   [0.0906391542748899, 0.0707225123213458, 0.0894257678938331, 0.0691401820113223, 0]])

# here we only plot the bottom left triangle, and the upper set into white
# set the xlim and ylim in 5
ax.set_xlim(0, 4)
ax.set_ylim(-1, 4)
# here we need switch the N2 and N1 
for i in range(5):
    for j in range(5):

        if i <=j:
            ax.text(j + 0.5, 5-i - 0.5, '', ha='center', va='center', color='white')
        else:
            color = plt.cm.viridis(matrix[i, j] / 0.1)
            val = matrix[i, j]

            if j == 0 and i > 1:
                box = plt.Rectangle((j+1, 5-i-1), 1, 1, color= color, ec='black')
                ax.add_patch(box)
                ax.text(j+1 + 0.5, 5-i - 0.5, f'{val:.3f}', ha='center', va='center', fontsize=16, color='white')

            elif j == 1:
                box = plt.Rectangle((j-1, 5-i-1), 1, 1, color= color, ec='black')
                ax.add_patch(box)
                ax.text(j-1 + 0.5, 5-i - 0.5, f'{val:.3f}', ha='center', va='center', fontsize=16, color='white')
            else:
                box = plt.Rectangle((j, 5-i-1), 1, 1, color= color, ec='black')
                ax.add_patch(box)
                ax.text(j + 0.5, 5-i - 0.5, f'{val:.3f}', ha='center', va='center', fontsize=16, color='white')

ax.text(-0.1, 0.5, 'South', ha='right', va='center', fontsize=14, rotation=90)
ax.text(-0.1, 1.5, 'North_4', ha='right', va='center', fontsize=14, rotation=90)
ax.text(-0.1, 2.5, 'North_3', ha='right', va='center', fontsize=14, rotation=90)
ax.text(-0.1, 3.5, 'North_2', ha='right', va='center', fontsize=14, rotation=90)

ax.text(0.5, -0.1, 'North_1', ha='center', va='top', fontsize=14)
ax.text(1.5, -0.1, 'North_2', ha='center', va='top', fontsize=14)
ax.text(2.5, -0.1, 'North_3', ha='center', va='top', fontsize=14)
ax.text(3.5, -0.1, 'North_4', ha='center', va='top', fontsize=14)
# off the axis
plt.axis('off')

# add the color bar legend, the max value is 0.1
cbar = plt.colorbar(plt.cm.ScalarMappable(cmap='viridis'), ax=ax, orientation='horizontal',pad=-0.08, aspect=50)
cbar.set_label('Pairwise Fst', fontsize=14)
# off the ticks
cbar.set_ticks([])
box = plt.Rectangle((0, -0.6), 0.02, 0.1, color='black', ec='black')
ax.add_patch(box)
plt.text(0, -0.5, '0', ha='center', va='bottom', fontsize=14, color='black')

box = plt.Rectangle((3.99, -0.6), 0.02, 0.1, color='black', ec='black')
ax.add_patch(box)
plt.text(4, -0.5, '0.10', ha='center', va='bottom', fontsize=12, color='black')

box = plt.Rectangle((2, -0.6), 0.02, 0.1, color='black', ec='black')
ax.add_patch(box)
plt.text(2, -0.5, '0.05', ha='center', va='bottom', fontsize=12, color='black')

plt.savefig('fst_matrix.png', dpi=300, bbox_inches='tight')

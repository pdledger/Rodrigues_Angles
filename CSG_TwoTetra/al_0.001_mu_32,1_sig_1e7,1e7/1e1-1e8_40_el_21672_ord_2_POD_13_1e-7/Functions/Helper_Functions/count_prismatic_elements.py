
def count_prismatic_elements(filename):
    """
    James Elgy - 2022
    Small function to count the number of prismatic elements in a mesh.
    This is done by evaluating the number of faces for each element in the vol file.
    :param filename: path for the .vol file.
    :return: number of prismatic elements and number of tetrahedral elements.
    """
    with open(filename, 'r') as f:
        stop = False
        line_number = 0
        while stop is False:
            line = f.readline()
            if line.rstrip() == 'volumeelements':
                stop = True
            line_number += 1

        max_elements = int(f.readline())

        stop = False
        line_number = 0
        while stop is False:
            line = f.readline()
            if len(line) > 2:
                if line[1:4] == ' 6 ':
                    stop = True
            line_number += 1
            if line[0:8] == '# surfid':
                stop = True
                line_number -= 3



        tet_elements = line_number
        prism_elements = max_elements - tet_elements

    return prism_elements, tet_elements

if __name__ == '__main__':
    filename = r'Results_08Nov/CSG_TwoTetra/al_0.001_mu_4,1_sig_1e7,1e7/1e1-1e8_40_el_21672_ord_3_POD_13_1e-7/Input_files/CSG_TwoTetra.vol'
    n_prisms, n_tets = count_prismatic_elements(filename)
    print(f' N Prisms = {n_prisms}, N Tets = {n_tets}')
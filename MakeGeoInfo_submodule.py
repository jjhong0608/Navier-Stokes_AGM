#!/usr/local/bin/python

def WriteRegion (fileID, regionName, conductivity, minX, minY, maxX, maxY, dx, dy, *args):
    if len(args) % 2 != 0: raise NameError ('argument number error')
    if len(args) < 1: raise NameError ("WriteRegion: There is no section")

    line = "%s\t%s\t" % ("REGION", regionName)
    fileID.write (line)

    line = "%23.16e\t%23.16e\t%23.16e\t%23.16e\t%23.16e\t%23.16e\t%23.16e\n" % (conductivity, minX, minY, maxX, maxY, dx, dy)
    fileID.write (line)

    for i in range (0, len (args), 2):
        fileID.write ("%s\t%s\n" % (args[i], args[i + 1]))

    fileID.write ("ENDREGION\n\n")

def WriteSectionLine (fileID, x1, y1, x2, y2, BC, BV):
    if BC == 'I':
        line = '%s\t%23.16e\t%23.16e\t%23.16e\t%23.16e\t%23.16e\t%23.16e\t%c\n' % ("LINE", x1, y1, x2, y2, y2 - y1, x1 - x2, BC)
    else:
        line = '%s\t%23.16e\t%23.16e\t%23.16e\t%23.16e\t%23.16e\t%23.16e\t%c\t%23.16e\n' % ("LINE", x1, y1, x2, y2, y2 - y1, x1 - x2, BC, BV)

    fileID.write (line)

def WriteLactangleSection (fileID, sectionName, minX, minY, maxX, maxY, BC, BV):
    if len (BC) != 4:
        raise NameError ('Boundary condition must have 4 elements')
    if len (BV) != 4:
        raise NameError ('Boundary value must have 4 elements')

    line = "%s\t%s\n" % ("SECTION", sectionName)
    fileID.write (line);

    WriteSectionLine (fileID, minX, minY, maxX, minY, BC[0], BV[0])
    WriteSectionLine (fileID, maxX, minY, maxX, maxY, BC[1], BV[1])
    WriteSectionLine (fileID, maxX, maxY, minX, maxY, BC[2], BV[2])
    WriteSectionLine (fileID, minX, maxY, minX, minY, BC[3], BV[3])

    fileID.write ("ENDSECTION\n\n")


def WriteOctagonSection (fileID, sectionName, minX, minY, maxX, maxY, BC, BV):
    if len (BC) != 8:
        raise NameError ('Boundary condition must have 8 elements')
    if len (BV) != 8:
        raise NameError ('Boundary value must have 8 elements')

    dx = (maxX - minX) / 3.0E0
    dy = (maxY - minY) / 3.0E0

    x = [minX, minX + dx, maxX - dx, maxX]
    y = [minY, minY + dy, maxY - dy, maxY]

    line = "%s\t%s\n" % ("SECTION", sectionName)
    fileID.write (line);

    WriteSectionLine (fileID, x[1], y[0], x[2], y[0], BC[0], BV[0])
    WriteSectionLine (fileID, x[2], y[0], x[3], y[1], BC[1], BV[1])
    WriteSectionLine (fileID, x[3], y[1], x[3], y[2], BC[2], BV[2])
    WriteSectionLine (fileID, x[3], y[2], x[2], y[3], BC[3], BV[3])
    WriteSectionLine (fileID, x[2], y[3], x[1], y[3], BC[4], BV[4])
    WriteSectionLine (fileID, x[1], y[3], x[0], y[2], BC[5], BV[5])
    WriteSectionLine (fileID, x[0], y[2], x[0], y[1], BC[6], BV[6])
    WriteSectionLine (fileID, x[0], y[1], x[1], y[0], BC[7], BV[7])

    fileID.write ("ENDSECTION\n\n")

def WriteCircleSection (fileID, sectionName, cX = 0.0E0, cY = 0.0E0, r = 1.0E0, BC = 'D', BV = 0.0E0, nLine = 101):
    import numpy as np

    theta = np.linspace (-np.pi, np.pi, nLine)
    x = r * np.cos (theta) + cX
    y = r * np.sin (theta) + cY

    line = "%s\t%s\n" % ("SECTION", sectionName)
    fileID.write (line);

    for i in range (len (theta) - 1):
        WriteSectionLine (fileID, x[i], y[i], x[i + 1], y[i + 1], BC, BV)

    fileID.write ("ENDSECTION\n\n")

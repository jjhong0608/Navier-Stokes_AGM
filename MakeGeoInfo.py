#!/usr/local/bin/python

def main():
    from MakeGeoInfo_submodule import WriteOctagonSection, WriteRegion, WriteCircleSection, WriteLactangleSection

    h = 0.1

    # translation_const = h * 0.5
    translation_const = 0.0E0
    with open("GeoInfo", 'w') as fileID:
        # WriteLactangleSection (fileID, "section0 0.0 0.0", 0.0E0, 0.0E0, 1.0E0, 1.0E0, ['D', 'D', 'D', 'D'], [0.0E0, 0.0E0, 0.0E0, 0.0E0])
        # WriteOctagonSection (fileID, "section0 0.0 0.0", -3.0E0, -3.0E0, 3.0E0, 3.0E0, ['N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'], [0.0E0, 0.0E0, 0.0E0, 0.0E0, 1.0E3, 0.0E0, 0.0E0, 0.0E0])
        # WriteOctagonSection (fileID, "section1 0.0 0.0", -1.5E0, -1.5E0, 1.5E0, 1.5E0, ['D', 'D', 'D', 'D', 'D', 'D', 'D', 'D'], [0.0E0, 0.0E0, 0.0E0, 0.0E0, 1.0E3, 0.0E0, 0.0E0, 0.0E0])
        WriteCircleSection(fileID, "section0 0.0 0.0", r=2.0E0, BC='N', BV=0, nLine=2001)
        WriteCircleSection(fileID, "section1 0.0 0.0", r=1.0E0, BC='D', BV=1e3, nLine=2001)

        # WriteRegion (fileID, "Region0", 1.0E0, 0.0E0, 0.0E0, 1.0E0, 1.0E0, h, h, "+", "section0")
        WriteRegion(fileID, "Region0", 5.0E-1, -2.0E0 - translation_const, -2.0E0 - translation_const,
                    2.0E0 + translation_const, 2.0E0 + translation_const, h, h, '+', "section0", '-', "section1")
        # WriteRegion (fileID, "Region0", 1.0E0, -3.0E0 - translation_const, -3.0E0 - translation_const / 3.14, 3.0E0 + translation_const, 3.0E0 + translation_const, h, h, '+', "section0", '-', "section1")
    # h = [0.1, 0.05, 0.025, 0.0125]
    # for i in range (4):
    #     with open ("GeoInfo" + str (i), 'w') as fileID:
    #         WriteOctagonSection (fileID, "section0 0.0 0.0", -3.0E0, -3.0E0, 3.0E0, 3.0E0, ['D', 'D', 'N', 'N', 'D', 'N', 'N', 'D'], [0.0E0, 0.0E0, 0.0E0, 0.0E0, 1.0E3, 0.0E0, 0.0E0, 0.0E0])
    #         # WriteOctagonSection (fileID, "section0 0.0 0.0", -3.0E0, -3.0E0, 3.0E0, 3.0E0, ['D', 'D', 'D', 'D', 'D', 'D', 'D', 'D'], [0.0E0, 0.0E0, 0.0E0, 0.0E0, 1.0E3, 0.0E0, 0.0E0, 0.0E0])
    #
    #         WriteRegion (fileID, "Region0", 1.0E0, -3.0E0 - translation_const, -3.0E0 - translation_const, 3.0E0 + translation_const, 3.0E0 + translation_const, h[i], h[i], '+', "section0")


if __name__ == "__main__":
    main()

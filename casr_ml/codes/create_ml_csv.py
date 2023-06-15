import sys
import os
import csv
import pandas as pd
import json
import operator
import re
import Bio
import math
import random
from Bio.Align import substitution_matrices
#@title Import Modules
import pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import ShuffleSplit



human_aa_dict = {1: 'M', 2: 'A', 3: 'F', 4: 'Y', 5: 'S', 6: 'C', 7: 'C', 8: 'W', 9: 'V', 10: 'L', 11: 'L', 12: 'A', 13: 'L', 14: 'T', 15: 'W', 16: 'H', 17: 'T', 18: 'S', 19: 'A', 20: 'Y', 21: 'G', 22: 'P', 23: 'D', 24: 'Q', 25: 'R', 26: 'A', 27: 'Q', 28: 'K', 29: 'K', 30: 'G', 31: 'D', 32: 'I', 33: 'I', 34: 'L', 35: 'G', 36: 'G', 37: 'L', 38: 'F', 39: 'P', 40: 'I', 41: 'H', 42: 'F', 43: 'G', 44: 'V', 45: 'A', 46: 'A', 47: 'K', 48: 'D', 49: 'Q', 50: 'D', 51: 'L', 52: 'K', 53: 'S', 54: 'R', 55: 'P', 56: 'E', 57: 'S', 58: 'V', 59: 'E', 60: 'C', 61: 'I', 62: 'R', 63: 'Y', 64: 'N', 65: 'F', 66: 'R', 67: 'G', 68: 'F', 69: 'R', 70: 'W', 71: 'L', 72: 'Q', 73: 'A', 74: 'M', 75: 'I', 76: 'F', 77: 'A', 78: 'I', 79: 'E', 80: 'E', 81: 'I', 82: 'N', 83: 'S', 84: 'S', 85: 'P', 86: 'A', 87: 'L', 88: 'L', 89: 'P', 90: 'N', 91: 'L', 92: 'T', 93: 'L', 94: 'G', 95: 'Y', 96: 'R', 97: 'I', 98: 'F', 99: 'D', 100: 'T', 101: 'C', 102: 'N', 103: 'T', 104: 'V', 105: 'S', 106: 'K', 107: 'A', 108: 'L', 109: 'E', 110: 'A', 111: 'T', 112: 'L', 113: 'S', 114: 'F', 115: 'V', 116: 'A', 117: 'Q', 118: 'N', 119: 'K', 120: 'I', 121: 'D', 122: 'S', 123: 'L', 124: 'N', 125: 'L', 126: 'D', 127: 'E', 128: 'F', 129: 'C', 130: 'N', 131: 'C', 132: 'S', 133: 'E', 134: 'H', 135: 'I', 136: 'P', 137: 'S', 138: 'T', 139: 'I', 140: 'A', 141: 'V', 142: 'V', 143: 'G', 144: 'A', 145: 'T', 146: 'G', 147: 'S', 148: 'G', 149: 'V', 150: 'S', 151: 'T', 152: 'A', 153: 'V', 154: 'A', 155: 'N', 156: 'L', 157: 'L', 158: 'G', 159: 'L', 160: 'F', 161: 'Y', 162: 'I', 163: 'P', 164: 'Q', 165: 'V', 166: 'S', 167: 'Y', 168: 'A', 169: 'S', 170: 'S', 171: 'S', 172: 'R', 173: 'L', 174: 'L', 175: 'S', 176: 'N', 177: 'K', 178: 'N', 179: 'Q', 180: 'F', 181: 'K', 182: 'S', 183: 'F', 184: 'L', 185: 'R', 186: 'T', 187: 'I', 188: 'P', 189: 'N', 190: 'D', 191: 'E', 192: 'H', 193: 'Q', 194: 'A', 195: 'T', 196: 'A', 197: 'M', 198: 'A', 199: 'D', 200: 'I', 201: 'I', 202: 'E', 203: 'Y', 204: 'F', 205: 'R', 206: 'W', 207: 'N', 208: 'W', 209: 'V', 210: 'G', 211: 'T', 212: 'I', 213: 'A', 214: 'A', 215: 'D', 216: 'D', 217: 'D', 218: 'Y', 219: 'G', 220: 'R', 221: 'P', 222: 'G', 223: 'I', 224: 'E', 225: 'K', 226: 'F', 227: 'R', 228: 'E', 229: 'E', 230: 'A', 231: 'E', 232: 'E', 233: 'R', 234: 'D', 235: 'I', 236: 'C', 237: 'I', 238: 'D', 239: 'F', 240: 'S', 241: 'E', 242: 'L', 243: 'I', 244: 'S', 245: 'Q', 246: 'Y', 247: 'S', 248: 'D', 249: 'E', 250: 'E', 251: 'E', 252: 'I', 253: 'Q', 254: 'H', 255: 'V', 256: 'V', 257: 'E', 258: 'V', 259: 'I', 260: 'Q', 261: 'N', 262: 'S', 263: 'T', 264: 'A', 265: 'K', 266: 'V', 267: 'I', 268: 'V', 269: 'V', 270: 'F', 271: 'S', 272: 'S', 273: 'G', 274: 'P', 275: 'D', 276: 'L', 277: 'E', 278: 'P', 279: 'L', 280: 'I', 281: 'K', 282: 'E', 283: 'I', 284: 'V', 285: 'R', 286: 'R', 287: 'N', 288: 'I', 289: 'T', 290: 'G', 291: 'K', 292: 'I', 293: 'W', 294: 'L', 295: 'A', 296: 'S', 297: 'E', 298: 'A', 299: 'W', 300: 'A', 301: 'S', 302: 'S', 303: 'S', 304: 'L', 305: 'I', 306: 'A', 307: 'M', 308: 'P', 309: 'Q', 310: 'Y', 311: 'F', 312: 'H', 313: 'V', 314: 'V', 315: 'G', 316: 'G', 317: 'T', 318: 'I', 319: 'G', 320: 'F', 321: 'A', 322: 'L', 323: 'K', 324: 'A', 325: 'G', 326: 'Q', 327: 'I', 328: 'P', 329: 'G', 330: 'F', 331: 'R', 332: 'E', 333: 'F', 334: 'L', 335: 'K', 336: 'K', 337: 'V', 338: 'H', 339: 'P', 340: 'R', 341: 'K', 342: 'S', 343: 'V', 344: 'H', 345: 'N', 346: 'G', 347: 'F', 348: 'A', 349: 'K', 350: 'E', 351: 'F', 352: 'W', 353: 'E', 354: 'E', 355: 'T', 356: 'F', 357: 'N', 358: 'C', 359: 'H', 360: 'L', 361: 'Q', 362: 'E', 363: 'G', 364: 'A', 365: 'K', 366: 'G', 367: 'P', 368: 'L', 369: 'P', 370: 'V', 371: 'D', 372: 'T', 373: 'F', 374: 'L', 375: 'R', 376: 'G', 377: 'H', 378: 'E', 379: 'E', 380: 'S', 381: 'G', 382: 'D', 383: 'R', 384: 'F', 385: 'S', 386: 'N', 387: 'S', 388: 'S', 389: 'T', 390: 'A', 391: 'F', 392: 'R', 393: 'P', 394: 'L', 395: 'C', 396: 'T', 397: 'G', 398: 'D', 399: 'E', 400: 'N', 401: 'I', 402: 'S', 403: 'S', 404: 'V', 405: 'E', 406: 'T', 407: 'P', 408: 'Y', 409: 'I', 410: 'D', 411: 'Y', 412: 'T', 413: 'H', 414: 'L', 415: 'R', 416: 'I', 417: 'S', 418: 'Y', 419: 'N', 420: 'V', 421: 'Y', 422: 'L', 423: 'A', 424: 'V', 425: 'Y', 426: 'S', 427: 'I', 428: 'A', 429: 'H', 430: 'A', 431: 'L', 432: 'Q', 433: 'D', 434: 'I', 435: 'Y', 436: 'T', 437: 'C', 438: 'L', 439: 'P', 440: 'G', 441: 'R', 442: 'G', 443: 'L', 444: 'F', 445: 'T', 446: 'N', 447: 'G', 448: 'S', 449: 'C', 450: 'A', 451: 'D', 452: 'I', 453: 'K', 454: 'K', 455: 'V', 456: 'E', 457: 'A', 458: 'W', 459: 'Q', 460: 'V', 461: 'L', 462: 'K', 463: 'H', 464: 'L', 465: 'R', 466: 'H', 467: 'L', 468: 'N', 469: 'F', 470: 'T', 471: 'N', 472: 'N', 473: 'M', 474: 'G', 475: 'E', 476: 'Q', 477: 'V', 478: 'T', 479: 'F', 480: 'D', 481: 'E', 482: 'C', 483: 'G', 484: 'D', 485: 'L', 486: 'V', 487: 'G', 488: 'N', 489: 'Y', 490: 'S', 491: 'I', 492: 'I', 493: 'N', 494: 'W', 495: 'H', 496: 'L', 497: 'S', 498: 'P', 499: 'E', 500: 'D', 501: 'G', 502: 'S', 503: 'I', 504: 'V', 505: 'F', 506: 'K', 507: 'E', 508: 'V', 509: 'G', 510: 'Y', 511: 'Y', 512: 'N', 513: 'V', 514: 'Y', 515: 'A', 516: 'K', 517: 'K', 518: 'G', 519: 'E', 520: 'R', 521: 'L', 522: 'F', 523: 'I', 524: 'N', 525: 'E', 526: 'E', 527: 'K', 528: 'I', 529: 'L', 530: 'W', 531: 'S', 532: 'G', 533: 'F', 534: 'S', 535: 'R', 536: 'E', 537: 'V', 538: 'P', 539: 'F', 540: 'S', 541: 'N', 542: 'C', 543: 'S', 544: 'R', 545: 'D', 546: 'C', 547: 'L', 548: 'A', 549: 'G', 550: 'T', 551: 'R', 552: 'K', 553: 'G', 554: 'I', 555: 'I', 556: 'E', 557: 'G', 558: 'E', 559: 'P', 560: 'T', 561: 'C', 562: 'C', 563: 'F', 564: 'E', 565: 'C', 566: 'V', 567: 'E', 568: 'C', 569: 'P', 570: 'D', 571: 'G', 572: 'E', 573: 'Y', 574: 'S', 575: 'D', 576: 'E', 577: 'T', 578: 'D', 579: 'A', 580: 'S', 581: 'A', 582: 'C', 583: 'N', 584: 'K', 585: 'C', 586: 'P', 587: 'D', 588: 'D', 589: 'F', 590: 'W', 591: 'S', 592: 'N', 593: 'E', 594: 'N', 595: 'H', 596: 'T', 597: 'S', 598: 'C', 599: 'I', 600: 'A', 601: 'K', 602: 'E', 603: 'I', 604: 'E', 605: 'F', 606: 'L', 607: 'S', 608: 'W', 609: 'T', 610: 'E', 611: 'P', 612: 'F', 613: 'G', 614: 'I', 615: 'A', 616: 'L', 617: 'T', 618: 'L', 619: 'F', 620: 'A', 621: 'V', 622: 'L', 623: 'G', 624: 'I', 625: 'F', 626: 'L', 627: 'T', 628: 'A', 629: 'F', 630: 'V', 631: 'L', 632: 'G', 633: 'V', 634: 'F', 635: 'I', 636: 'K', 637: 'F', 638: 'R', 639: 'N', 640: 'T', 641: 'P', 642: 'I', 643: 'V', 644: 'K', 645: 'A', 646: 'T', 647: 'N', 648: 'R', 649: 'E', 650: 'L', 651: 'S', 652: 'Y', 653: 'L', 654: 'L', 655: 'L', 656: 'F', 657: 'S', 658: 'L', 659: 'L', 660: 'C', 661: 'C', 662: 'F', 663: 'S', 664: 'S', 665: 'S', 666: 'L', 667: 'F', 668: 'F', 669: 'I', 670: 'G', 671: 'E', 672: 'P', 673: 'Q', 674: 'D', 675: 'W', 676: 'T', 677: 'C', 678: 'R', 679: 'L', 680: 'R', 681: 'Q', 682: 'P', 683: 'A', 684: 'F', 685: 'G', 686: 'I', 687: 'S', 688: 'F', 689: 'V', 690: 'L', 691: 'C', 692: 'I', 693: 'S', 694: 'C', 695: 'I', 696: 'L', 697: 'V', 698: 'K', 699: 'T', 700: 'N', 701: 'R', 702: 'V', 703: 'L', 704: 'L', 705: 'V', 706: 'F', 707: 'E', 708: 'A', 709: 'K', 710: 'I', 711: 'P', 712: 'T', 713: 'S', 714: 'F', 715: 'H', 716: 'R', 717: 'K', 718: 'W', 719: 'W', 720: 'G', 721: 'L', 722: 'N', 723: 'L', 724: 'Q', 725: 'F', 726: 'L', 727: 'L', 728: 'V', 729: 'F', 730: 'L', 731: 'C', 732: 'T', 733: 'F', 734: 'M', 735: 'Q', 736: 'I', 737: 'V', 738: 'I', 739: 'C', 740: 'V', 741: 'I', 742: 'W', 743: 'L', 744: 'Y', 745: 'T', 746: 'A', 747: 'P', 748: 'P', 749: 'S', 750: 'S', 751: 'Y', 752: 'R', 753: 'N', 754: 'Q', 755: 'E', 756: 'L', 757: 'E', 758: 'D', 759: 'E', 760: 'I', 761: 'I', 762: 'F', 763: 'I', 764: 'T', 765: 'C', 766: 'H', 767: 'E', 768: 'G', 769: 'S', 770: 'L', 771: 'M', 772: 'A', 773: 'L', 774: 'G', 775: 'F', 776: 'L', 777: 'I', 778: 'G', 779: 'Y', 780: 'T', 781: 'C', 782: 'L', 783: 'L', 784: 'A', 785: 'A', 786: 'I', 787: 'C', 788: 'F', 789: 'F', 790: 'F', 791: 'A', 792: 'F', 793: 'K', 794: 'S', 795: 'R', 796: 'K', 797: 'L', 798: 'P', 799: 'E', 800: 'N', 801: 'F', 802: 'N', 803: 'E', 804: 'A', 805: 'K', 806: 'F', 807: 'I', 808: 'T', 809: 'F', 810: 'S', 811: 'M', 812: 'L', 813: 'I', 814: 'F', 815: 'F', 816: 'I', 817: 'V', 818: 'W', 819: 'I', 820: 'S', 821: 'F', 822: 'I', 823: 'P', 824: 'A', 825: 'Y', 826: 'A', 827: 'S', 828: 'T', 829: 'Y', 830: 'G', 831: 'K', 832: 'F', 833: 'V', 834: 'S', 835: 'A', 836: 'V', 837: 'E', 838: 'V', 839: 'I', 840: 'A', 841: 'I', 842: 'L', 843: 'A', 844: 'A', 845: 'S', 846: 'F', 847: 'G', 848: 'L', 849: 'L', 850: 'A', 851: 'C', 852: 'I', 853: 'F', 854: 'F', 855: 'N', 856: 'K', 857: 'I', 858: 'Y', 859: 'I', 860: 'I', 861: 'L', 862: 'F', 863: 'K', 864: 'P', 865: 'S', 866: 'R', 867: 'N', 868: 'T', 869: 'I', 870: 'E', 871: 'E', 872: 'V', 873: 'R', 874: 'C', 875: 'S', 876: 'T', 877: 'A', 878: 'A', 879: 'H', 880: 'A', 881: 'F', 882: 'K', 883: 'V', 884: 'A', 885: 'A', 886: 'R', 887: 'A', 888: 'T', 889: 'L', 890: 'R', 891: 'R', 892: 'S', 893: 'N', 894: 'V', 895: 'S', 896: 'R', 897: 'K', 898: 'R', 899: 'S', 900: 'S', 901: 'S', 902: 'L', 903: 'G', 904: 'G', 905: 'S', 906: 'T', 907: 'G', 908: 'S', 909: 'T', 910: 'P', 911: 'S', 912: 'S', 913: 'S', 914: 'I', 915: 'S', 916: 'S', 917: 'K', 918: 'S', 919: 'N', 920: 'S', 921: 'E', 922: 'D', 923: 'P', 924: 'F', 925: 'P', 926: 'Q', 927: 'P', 928: 'E', 929: 'R', 930: 'Q', 931: 'K', 932: 'Q', 933: 'Q', 934: 'Q', 935: 'P', 936: 'L', 937: 'A', 938: 'L', 939: 'T', 940: 'Q', 941: 'Q', 942: 'E', 943: 'Q', 944: 'Q', 945: 'Q', 946: 'Q', 947: 'P', 948: 'L', 949: 'T', 950: 'L', 951: 'P', 952: 'Q', 953: 'Q', 954: 'Q', 955: 'R', 956: 'S', 957: 'Q', 958: 'Q', 959: 'Q', 960: 'P', 961: 'R', 962: 'C', 963: 'K', 964: 'Q', 965: 'K', 966: 'V', 967: 'I', 968: 'F', 969: 'G', 970: 'S', 971: 'G', 972: 'T', 973: 'V', 974: 'T', 975: 'F', 976: 'S', 977: 'L', 978: 'S', 979: 'F', 980: 'D', 981: 'E', 982: 'P', 983: 'Q', 984: 'K', 985: 'N', 986: 'A', 987: 'M', 988: 'A', 989: 'H', 990: 'R', 991: 'N', 992: 'S', 993: 'T', 994: 'H', 995: 'Q', 996: 'N', 997: 'S', 998: 'L', 999: 'E', 1000: 'A', 1001: 'Q', 1002: 'K', 1003: 'S', 1004: 'S', 1005: 'D', 1006: 'T', 1007: 'L', 1008: 'T', 1009: 'R', 1010: 'H', 1011: 'E', 1012: 'P', 1013: 'L', 1014: 'L', 1015: 'P', 1016: 'L', 1017: 'Q', 1018: 'C', 1019: 'G', 1020: 'E', 1021: 'T', 1022: 'D', 1023: 'L', 1024: 'D', 1025: 'L', 1026: 'T', 1027: 'V', 1028: 'Q', 1029: 'E', 1030: 'T', 1031: 'G', 1032: 'L', 1033: 'Q', 1034: 'G', 1035: 'P', 1036: 'V', 1037: 'G', 1038: 'G', 1039: 'D', 1040: 'Q', 1041: 'R', 1042: 'P', 1043: 'E', 1044: 'V', 1045: 'E', 1046: 'D', 1047: 'P', 1048: 'E', 1049: 'E', 1050: 'L', 1051: 'S', 1052: 'P', 1053: 'A', 1054: 'L', 1055: 'V', 1056: 'V', 1057: 'S', 1058: 'S', 1059: 'S', 1060: 'Q', 1061: 'S', 1062: 'F', 1063: 'V', 1064: 'I', 1065: 'S', 1066: 'G', 1067: 'G', 1068: 'G', 1069: 'S', 1070: 'T', 1071: 'V', 1072: 'T', 1073: 'E', 1074: 'N', 1075: 'V', 1076: 'V', 1077: 'N', 1078: 'S'}

tm_tendency_dict = {"A":0.380 , "R":-2.570,"N":-1.620,"D":-3.270 , "C":-0.300 ,  "Q":-1.840 , "E": -2.900 , "G":-0.190  ,
"H":-1.440, "I":1.970,  "L":1.820,  "K":-3.460 , "M":1.400 ,  "F":1.980,  "P":-1.440 ,  "S":-0.530  ,  "T":-0.320 ,  
"W":1.530   ,  "Y":0.490  ,  "V":1.460}


molecular_weight_dict = {"A": 89.000, "R":174.000,"N":132.000,"D": 133.000, "C": 121.000,  "Q": 146.000, "E": 147.000, "G": 75.000,
"H":155.000, "I":131.000,  "L":131.000,  "K": 146.000, "M": 149.000,  "F":165.000,  "P": 115.000,  "S": 105.000,  "T": 119.000,  
"W": 204.000,  "Y": 181.000,  "V": 117.000}

bulkiness_dict = {"A": 11.500, "R":14.280,"N":12.820,"D": 11.680, "C": 13.460,  "Q": 14.450, "E": 13.570, "G": 3.400,
"H":13.690, "I":21.400,  "L":21.400,  "K": 15.710, "M": 16.250,  "F":19.800,  "P": 17.430,  "S": 9.470,  "T": 15.770,  
"W": 21.670,  "Y": 18.030,  "V": 21.570}



atomic_weight_ratio_dict = {"A":0.000  ,"R":0.650, "N":1.330 ,"D":1.380   ,"C":2.750  ,"Q":0.890 ,"E":0.920  ,"G":0.740 ,"H":0.580,  
"I":0.000  ,  "L":0.000,"K":0.330, "M":0.000 ,"F":0.000 ,"P":0.390,"S":1.420   ,"T":0.710  ,"W":0.130 ,  "Y":0.200 , "V":0.000}

doolittle_hydropathicity_dict =  {"A":1.800 ,"R":-4.500, "N":-3.500 ,"D":-3.500 ,"C":2.500  ,"Q":-3.500,"E":-3.500 ,"G":-0.400 ,"H":-3.200,  
"I":4.500  ,  "L":3.800,"K":-3.900, "M":1.900 ,"F":2.800 ,"P":-1.600,"S":-0.800  ,"T":-0.700  ,"W":-0.900 ,  "Y":-1.300,  "V":4.200}


zimmerman_polarity_dict = {"A":0.000,"R":52.000, "N": 3.380,"D": 49.700,"C":1.480,"Q":3.530,"E":49.900,"G":0.000,"H":51.600,  
"I":0.130,  "L":0.130,"K":49.500, "M":1.430,"F":0.350,"P":1.580,"S":1.670,"T":1.660,"W":2.100,  "Y":1.610,  "V":0.130}

average_flexibility_dict = {"A":0.360,"R":0.530,"N":0.460,"D":0.510,"C":0.350,"Q":0.490,"E":0.500, "G":0.540,"H":0.320,
"I":0.460, "L":0.370,  "K":0.470,  "M":0.300,  "F":0.310,"P":0.510,   "S":0.510,   "T":0.440,  "W":0.310,  "Y":0.420 ,"V":0.390}

dayhoff_dict = {"A":100.000,"R":65.000,  "N":134.000,"D":106.000,"C":20.000,"Q":93.000, "E":102.000, "G":49.000,
"H":66.000,  "I":96.000,  "L":40.000, "K":56.000,  "M":94.000, "F":41.000,  "P":56.000,"S":120.000,  "T":97.000,  
"W":18.000, "Y":41.000, "V":74.000}  

avg_buried_area_dict ={"A":86.600,"R":162.200,"N":103.300,"D":97.800, "C":132.300, "Q":119.200,"E":113.900,"G": 62.900,  
"H": 155.800, "I":158.000,"L":164.100,"K":115.500,"M":172.900,"F":194.100,"P":92.900,"S":85.600,"T":106.500,"W":224.600,  
"Y": 177.700,"V": 141.000} 


aa_snc_dict = {"Gly":"G", "Ala":"A","Val":"V", "Leu":"L","Ile":"I","Met":"M","Phe":"F","Trp":"W","Pro":"P",
"Ser":"S","Thr":"T","Cys":"C","Tyr":"Y","Asn":"N","Gln":"Q",
"Asp":"D","Glu":"E",
"Lys":"K","Arg":"R","His":"H"}

disease_cause_dict = {"FHH":"loss-of-function", "NSHPT":"loss-of-function", "ADH":"gain-of-function","Hyperparathyroidism":"loss-of-function",
"NSHPT, FHH":"loss-of-function", "Hypoparathyroidism":"gain-of-function"}

mutations_dict = {"p.Met1Arg":"NSHPT",
"p.Leu11Ser":"FHH",
"p.Leu13Pro":"FHH",
"p.Gly27Pro":"FHH",
"p.Gln27Glu":"ADH",
"p.Gln27Arg":"FHH",
"p.Lys29Glu":"ADH",
"p.Pro39Ala":"FHH",
"p.Ile40Phe":"FHH",
"p.Lys47Asn":"ADH",
"p.Ser53Pro":"FHH",
"p.Pro55Leu":"FHH",
"p.Arg62Met":"NSHPT, FHH",
"p.Arg66His":"NSHPT, FHH",
"p.Arg66Cys":"NSHPT, FHH",
"p.Ile81Met":"FHH",
"p.Thr100Ile":"Hyperparathyroidism",
"p.Ala110Thr":"FHH",
"p.Ala116Thr":"ADH",
"p.Asn118Lys":"ADH",
"p.Ser122Cys":"ADH",
"p.Leu123Ser":"ADH",
"p.Asn124Lys":"ADH",
"p.Leu125Pro":"ADH",
"p.Asp126Val":"ADH",
"p.Phe128Leu":"ADH",
"p.Cys129Ser":"ADH",
"p.Cys131Trp":"ADH",
"p.Cys131Ser":"ADH",
"p.Pro136Leu":"ADH",
"p.Ser137Pro":"FHH",
"p.Thr138Met":"FHH",
"p.Gly143Glu":"FHH",
"p.Thr145Ile":"FHH",
"p.Thr151Arg":"ADH",
"p.Thr151Met":"Hypoparathyroidism",
"p.Leu159Phe":"FHH",
"p.Leu159Pro":"Hyperparathyroidism",
"p.Arg172Gly":"FHH",
"p.Leu173Pro":"FHH",
"p.Leu173Phe":"ADH",
"p.Asn178Tyr":"ADH",
"p.Asn178Asp":"FHH",
"p.Phe180Cys":"FHH",
"p.Arg185Gln":"NSHPT, FHH",
"p.Asp190Gly":"FHH",
"p.Glu191Lys":"ADH",
"p.Ile212Thr":"FHH",
"p.Asp215Gly":"FHH",
"p.Tyr218Cys":"FHH",
"p.Arg220Trp":"NSHPT, FHH",
"p.Pro221Leu":"ADH",
"p.Pro221Gln":"FHH",
"p.Arg227Gly":"FHH",
"p.Arg227Gln":"Hyperparathyroidism",
"p.Arg227Leu":"NSHPT",
"p.Glu228Gln":"ADH",
"p.Glu241Lys":"ADH",
"p.Gln245Arg":"ADH",
"p.Thr263Met":"FHH",
"p.Ser272Ile":"FHH",
"p.Ile283Thr":"Hyperparathyroidism",
"p.Glu297Lys":"NSHPT, FHH",
"p.Glu297Asp":"ADH",
"p.Pro339Thr":"Hyperparathyroidism",
"p.Phe351Val":"FHH",
"p.Cys395Arg":"FHH",
"p.Asp410Glu":"ADH",
"p.Ser417Cys":"FHH",
"p.Gln459Arg":"FHH",
"p.Arg465Trp":"FHH",
"p.Arg465Gln":"FHH",
"p.Glu481Lys":"ADH",
"p.Trp530Gly":"FHH",
"p.Arg544Gln":"ADH",
"p.Gly549Arg":"FHH",
"p.Thr550Ile":"FHH",
"p.Arg551Lys":"NSHPT",
"p.Gly553Arg":"FHH",
"p.Gly557Glu":"FHH",
"p.Cys568Tyr":"FHH",
"p.Pro569His":"ADH",
"p.Gly571Trp":"FHH",
"p.Gly571Val":"FHH",
"p.Cys582Arg":"FHH",
"p.Cys582Tyr":"NSHPT",
"p.His595Tyr":"FHH",
"p.Glu604Lys":"ADH",
"p.Phe612Ser":"ADH",
"p.Leu616Val":"ADH",
"p.Arg638Leu":"FHH",
"p.Leu650Pro":"Hyperparathyroidism",
"p.Ser657Tyr":"FHH",
"p.Gly670Glu":"NSHPT",
"p.Gly670Arg":"FHH",
"p.Arg680Cys":"NSHPT, FHH",
"p.Arg680Gly":"ADH",
"p.Arg680His":"FHH",
"p.Gln681His":"ADH",
"p.Gln681Arg":"ADH",
"p.Val689Met":"Hyperparathyroidism",
"p.Val697Met":"FHH",
"p.Glu707Val":"FHH",
"p.Leu727Gln":"ADH",
"p.Met734Arg":"FHH",
"p.Pro748His":"FHH",
"p.Pro748Arg":"FHH",
"p.Cys765Trp":"FHH",
"p.Glu767Lys":"ADH",
"p.Glu767Gln":"ADH",
"p.Leu773Arg":"Hypoparathyroidism",
"p.Gly774Ser":"FHH",
"p.Ile777Leu":"FHH",
"p.Gly778Asp":"FHH",
"p.Phe788Leu":"ADH",
"p.Phe788Cys":"Hypoparathyroidism",
"p.Arg795Trp":"FHH",
"p.Arg795Gln":"FHH",
"p.Arg795Pro":"FHH",
"p.Pro798Thr":"FHH",
"p.Glu799Lys":"Hypoparathyroidism",
"p.Asn802Ile":"ADH",
"p.Asn802Ser":"FHH",
"p.Ala804Asp":"FHH",
"p.Phe806Ser":"ADH",
"p.Val817Ile":"FHH",
"p.Ser820Phe":"ADH",
"p.Ser820Ala":"FHH",
"p.Phe821Leu":"Hypoparathyroidism",
"p.Ala824Pro":"ADH",
"p.Ala824Ser":"Hypoparathyroidism",
"p.Thr828Asn":"ADH",
"p.Thr828Ile":"FHH",
"p.Gly830Ser":"ADH",
"p.Phe832Ser":"ADH",
"p.Ala835Thr":"ADH",
"p.Val836Leu":"Hypoparathyroidism",
"p.Glu837Asp":"ADH",
"p.Ile839Thr":"ADH",
"p.Ala843Glu":"ADH",
"p.Ala844Thr":"ADH",
"p.Leu849Pro":"FHH",
"p.Cys851Phe":"FHH",
"p.Phe881Leu":"FHH",
"p.Thr888Met":"ADH",
"p.Gln926Arg":"FHH",
"p.Thr972Met":"FHH",
"p.Asp1005Asn":"FHH",
"p.Gly1019Arg":"FHH",
"p.Ala19Pro":"FHH",
"p.Gly21Arg":"FHH",
"p.Ala26Gly":"FHH",
"p.Ala26Ser":"FHH",
"p.Asp31Gly":"FHH",
"p.Leu37Pro":"FHH",
"p.Phe42Ser":"FHH",
"p.Leu51Arg":"FHH",
"p.Cys60Phe":"NSHPT, FHH",
"p.Cys60Arg":"FHH",
"p.Tyr63Cys":"FHH",
"p.Gly67Val":"FHH",
"p.Arg69Cys":"FHH",
"p.Arg69His":"NSHPT, FHH",
"p.Met74Val":"FHH",
"p.Met74Thr":"FHH",
"p.Ala77Val":"FHH",
"p.Glu80Asp":"FHH",
"p.Ile81Thr":"FHH",
"p.Ile81Lys":"NSHPT",
"p.Asn82Lys":"FHH",
"p.Asn90Thr":"FHH",
"p.Cys101Arg":"FHH",
"p.Cys101Trp":"FHH",
"p.Val104Ile":"ADH",
"p.Val104Pro":"FHH",
"p.Asn116Pro":"ADH",
"p.Lys119Ile":"ADH",
"p.Leu123Met":"ADH",
"p.Leu125Phe":"Hypoparathyroidism",
"p.Glu127Lys":"Hypoparathyroidism",
"p.Glu127Gly":"ADH",
"p.Glu127Ala":"ADH",
"p.Phe128Ala":"ADH",
"p.Cys129Phe":"ADH",
"p.Cys129Arg":"ADH",
"p.Cys129Tyr":"ADH",
"p.Cys131Tyr":"ADH",
"p.Cys131Phe":"Hypoparathyroidism",
"p.Glu133Val":"ADH",
"p.Ala140Val":"FHH",
"p.Gly143Arg":"FHH",
"p.Gly146Asp":"FHH",
"p.Val149Leu":"FHH",
"p.Ser150Cys":"FHH",
"p.Gly158Glu":"FHH",
"p.Gly158Arg":"FHH",
"p.Leu159Arg":"FHH",
"p.Tyr161Cys":"Hyperparathyroidism",
"p.Ile162Phe":"FHH",
"p.Pro163Arg":"FHH",
"p.Gln164Lys":"FHH",
"p.Val165Asp":"FHH",
"p.Ser166Gly":"FHH",
"p.Ser171Arg":"FHH",
"p.Ser171Asn":"FHH",
"p.Leu174Arg":"FHH",
"p.Ser175Ile":"FHH",
"p.Asn178Thr":"FHH",
"p.Phe183Leu":"FHH",
"p.Arg185Gly":"FHH",
"p.Asp190Glu":"FHH",
"p.Arg205Cys":"FHH",
"p.Trp208Ser":"FHH",
"p.Trp208Cys":"FHH",
"p.Ile212Ser":"NSHPT",
"p.Ala213Glu":"FHH",
"p.Tyr218His":"FHH",
"p.Tyr218Ser":"FHH",
"p.Arg220Pro":"FHH",
"p.Arg220Gln":"FHH",
"p.Pro221Ser":"FHH",
"p.Gly222Glu":"FHH",
"p.Lys225Thr":"FHH",
"p.Glu228Lys":"ADH",
"p.Glu228Gly":"ADH",
"p.Glu231Lys":"FHH",
"p.Cys236Gly":"FHH",
"p.Gln253His":"FHH",
"p.Ser271Phe":"FHH",
"p.Ile288Val":"FHH",
"p.Ser296Asn":"FHH",
"p.Ala298Thr":"Hyperparathyroidism",
"p.Ser302Phe":"FHH",
"p.Pro339Ser":"FHH",
"p.Asn345Ile":"FHH",
"p.Phe351Ile":"FHH",
"p.Cys358Ser":"FHH",
"p.Gly397Arg":"FHH",
"p.Tyr418Cys":"FHH",
"p.Asn419Ser":"ADH",
"p.Ala423Lys":"FHH",
"p.Ala430Val":"FHH",
"p.Gln432Arg":"FHH",
"p.Cys437Arg":"FHH",
"p.Ala457Val":"FHH",
"p.Leu461Pro":"FHH",
"p.His463Gln":"FHH",
"p.Phe469Cys":"FHH",
"p.Asn493Lys":"FHH",
"p.Val504Met":"FHH",
"p.Gly509Glu":"FHH",
"p.Gly509Arg":"NSHPT, FHH",
"p.Gly509Ala":"FHH",
"p.Tyr514Cys":"Hypoparathyroidism",
"p.Glu536Asp":"FHH",
"p.Cys542Tyr":"FHH",
"p.Cys542Arg":"FHH",
"p.Cys546Ser":"FHH",
"p.Cys546Gly":"FHH",
"p.Arg551Gly":"FHH",
"p.Ile555Thr":"NSHPT, FHH",
"p.Ile555Val":"FHH",
"p.Glu556Lys":"ADH",
"p.Cys561Gly":"FHH",
"p.Cys562Tyr":"FHH",
"p.Cys565Gly":"FHH",
"p.Asp570Gly":"FHH",
"p.Tyr573Cys":"FHH",
"p.Cys582Phe":"FHH",
"p.Cys585Phe":"FHH",
"p.Phe589Leu":"ADH",
"p.Ser591Cys":"NSHPT, FHH",
"p.Lys601Asn":"FHH",
"p.Glu610Gly":"ADH",
"p.Gly613Arg":"FHH",
"p.Gly623Asp":"FHH",
"p.Phe634Ser":"FHH",
"p.Thr640Ile":"FHH",
"p.Ala645Asp":"FHH",
"p.Cys661Tyr":"FHH",
"p.Cys661Trp":"FHH",
"p.Phe668Leu":"FHH",
"p.Cys677Tyr":"FHH",
"p.Arg678His":"FHH",
"p.Leu690Arg":"FHH",
"p.Ser693Leu":"FHH",
"p.Arg701Pro":"FHH",
"p.Leu704Pro":"Hyperparathyroidism",
"p.Arg716His":"FHH",
"P.Arg716Cys":"FHH",
"p.Asn722Ser":"ADH",
"p.Val728Phe":"FHH",
"p.Met734Thr":"Hypoparathyroidism",
"p.Gln735Pro":"ADH",
"p.Trp742Arg":"FHH",
"p.Leu743Pro":"FHH",
"p.Ala746Val":"FHH",
"p.Pro747Leu":"FHH",
"p.Pro748Leu":"FHH",
"p.Arg752Cys":"FHH",
"p.Glu757Lys":"ADH",
"p.Ile760Asn":"FHH",
"p.Cys765Arg":"FHH",
"p.Gly768Val":"NSHPT",
"p.Gly778Arg":"FHH",
"p.Cys781Arg":"FHH",
"p.Leu782Arg":"FHH",
"p.Ala791Pro":"FHH",
"p.Pro798Leu":"FHH",
"p.Asn800Thr":"ADH",
"p.Phe801Leu":"FHH",
"p.Ile807Asn":"FHH",
"p.Phe809Ile":"FHH",
"p.Phe809Leu":"FHH",
"p.Met811Val":"ADH",
"p.Phe816Thr":"FHH",
"p.Phe816Val":"FHH",
"p.Trp818Leu":"ADH",
"p.Tyr825Pro":"ADH",
"p.Ser827Arg":"FHH",
"p.Tyr829Cys":"ADH",
"p.Lys831Thr":"ADH",
"p.Phe832Leu":"ADH",
"p.Ala844Pro":"ADH",
"p.Ser845Asn":"ADH",
"p.Ile860Ser":"FHH",
"p.Leu861Pro":"FHH",
"p.Asn867Ile":"FHH",
"p.Ile869Met":"FHH",
"p.Arg873His":"FHH",
"p.Arg873Pro":"FHH",
"p.Val883Met":"ADH",
"p.Ala885Asp":"FHH",
"p.Arg886Trp":"NSHPT, FHH",
"p.Arg886Pro":"Hyperparathyroidism",
"p.Ala887Asp":"FHH",
"p.Asn1074Asp":"FHH"}

neutral_subs_msa_dict = {2: {'T': 0.3024193548387097}, 3: {'L': 0.42338709677419356}, 5: {'L': 0.21138211382113822}, 7: {'Y': 0.22764227642276422}, 8: {'L': 0.7479674796747967}, 9: {'I': 0.5853658536585366}, 12: {'G': 0.20647773279352227, 'L': 0.20242914979757085}, 13: {'F': 0.5506072874493927}, 16: {'V': 0.20242914979757085, 'N': 0.2550607287449393}, 18: {'A': 0.22672064777327935}, 23: {'N': 0.3024193548387097}, 29: {'T': 0.208}, 33: {'L': 0.21052631578947367}, 46: {'S': 0.2289156626506024}, 52: {'A': 0.21686746987951808}, 53: {'A': 0.21370967741935484}, 58: {'T': 0.2289156626506024}, 61: {'V': 0.22983870967741934}, 83: {'N': 0.2749003984063745}, 85: {'S': 0.21031746031746032}, 86: {'T': 0.4246031746031746}, 91: {'M': 0.6573705179282868, 'I': 0.23904382470119523}, 132: {'T': 0.21912350597609562}, 133: {'D': 0.2222222222222222}, 149: {'I': 0.5577689243027888}, 165: {'I': 0.21825396825396826}, 180: {'Y': 0.2698412698412698}, 184: {'M': 0.2222222222222222}, 189: {'T': 0.20238095238095238}, 205: {'Q': 0.21031746031746032}, 210: {'I': 0.2222222222222222}, 211: {'A': 0.2222222222222222}, 214: {'S': 0.2222222222222222}, 227: {'E': 0.2222222222222222}, 230: {'M': 0.20238095238095238}, 239: {'L': 0.2222222222222222}, 240: {'N': 0.20238095238095238}, 254: {'Q': 0.6111111111111112}, 255: {'L': 0.21825396825396826}, 257: {'D': 0.2261904761904762}, 258: {'R': 0.21825396825396826}, 260: {'E': 0.20634920634920634}, 265: {'R': 0.25793650793650796}, 271: {'A': 0.21825396825396826}, 283: {'M': 0.21031746031746032}, 290: {'D': 0.21825396825396826}, 291: {'R': 0.6349206349206349}, 307: {'K': 0.22134387351778656}, 309: {'E': 0.9365079365079365}, 310: {'F': 0.3333333333333333}, 311: {'L': 0.21031746031746032}, 312: {'D': 0.21428571428571427, 'R': 0.2261904761904762}, 314: {'I': 0.28286852589641437}, 316: {'S': 0.24206349206349206}, 335: {'Q': 0.8650793650793651}, 340: {'K': 0.5140562248995983}, 343: {'S': 0.21428571428571427, 'T': 0.23015873015873015}, 344: {'N': 0.28286852589641437}, 346: {'E': 0.21031746031746032}, 348: {'V': 0.23412698412698413}, 349: {'R': 0.21825396825396826}, 359: {'Y': 0.5198412698412699}, 361: {'P': 0.30039525691699603}, 362: {'D': 0.4624505928853755}, 363: {'S': 0.22924901185770752}, 364: {'S': 0.26877470355731226}, 366: {'N': 0.3384615384615385}, 367: {'S': 0.3384615384615385}, 368: {'P': 0.35751295336787564}, 369: {'A': 0.32642487046632124}, 370: {'S': 0.3492063492063492, 'M': 0.2328042328042328}, 371: {'T': 0.2}, 372: {'S': 0.38144329896907214}, 374: {'H': 0.36082474226804123}, 375: {'K': 0.324}, 380: {'G': 0.5983606557377049}, 381: {'L': 0.27835051546391754}, 382: {'G': 0.6649746192893401}, 384: {'A': 0.2849740932642487, 'I': 0.22797927461139897}, 385: {'G': 0.37823834196891193}, 387: {'G': 0.456}, 388: {'T': 0.30158730158730157}, 389: {'A': 0.21115537848605578}, 390: {'S': 0.20967741935483872}, 394: {'P': 0.30278884462151395}, 398: {'E': 0.32142857142857145}, 400: {'D': 0.2222222222222222}, 402: {'T': 0.43253968253968256}, 409: {'M': 0.6428571428571429, 'L': 0.2619047619047619}, 422: {'V': 0.21031746031746032}, 438: {'T': 0.48412698412698413}, 441: {'K': 0.3888888888888889}, 445: {'A': 0.2222222222222222}, 472: {'S': 0.23137254901960785}, 478: {'D': 0.3425196850393701}, 482: {'F': 0.2755905511811024}, 490: {'T': 0.2283464566929134}, 496: {'R': 0.21653543307086615}, 503: {'V': 0.49606299212598426}, 506: {'E': 0.44881889763779526}, 510: {'H': 0.33070866141732286}, 524: {'D': 0.2047244094488189}, 526: {'N': 0.2874015748031496}, 531: {'N': 0.23137254901960785}, 535: {'K': 0.2795275590551181}, 544: {'E': 0.2178988326848249}, 548: {'P': 0.5348837209302325}, 566: {'T': 0.21621621621621623}, 567: {'D': 0.2548262548262548}, 569: {'S': 0.2364341085271318}, 576: {'H': 0.2140077821011673}, 581: {'V': 0.2099236641221374}, 583: {'D': 0.6297709923664122}, 587: {'E': 0.2813688212927757, 'N': 0.22053231939163498}, 588: {'N': 0.21292775665399238}, 589: {'S': 0.23574144486692014, 'Y': 0.25475285171102663}, 599: {'F': 0.21292775665399238}, 600: {'P': 0.3041825095057034, 'L': 0.20912547528517111}, 602: {'Q': 0.34220532319391633}, 624: {'V': 0.20152091254752852}, 628: {'S': 0.3574144486692015}, 635: {'T': 0.3041825095057034}, 659: {'I': 0.23954372623574144}, 674: {'N': 0.23954372623574144}, 714: {'L': 0.5381679389312977}, 734: {'V': 0.5610687022900763}, 736: {'V': 0.21374045801526717}, 737: {'M': 0.21374045801526717}, 741: {'V': 0.22433460076045628}, 745: {'N': 0.23574144486692014}, 754: {'H': 0.8631178707224335}, 755: {'D': 0.24427480916030533}, 756: {'I': 0.21374045801526717}, 766: {'N': 0.25475285171102663}, 802: {'T': 0.21292775665399238}, 826: {'F': 0.21374045801526717}, 844: {'S': 0.2480916030534351}, 857: {'V': 0.8778625954198473}, 886: {'K': 0.21621621621621623}, 900: {'N': 0.4163424124513619}, 906: {'S': 0.20463320463320464}, 918: {'T': 0.2140077821011673}, 920: {'H': 0.308}, 926: {'L': 0.2890625}, 928: {'A': 0.359375}, 940: {'S': 0.46502057613168724}, 941: {'A': 0.29508196721311475}, 951: {'Q': 0.30303030303030304}, 954: {'R': 0.23863636363636365}, 955: {'Q': 0.5454545454545454}, 956: {'R': 0.25}, 960: {'R': 0.3891402714932127}, 961: {'G': 0.3738738738738739}, 964: {'P': 0.22178988326848248}, 965: {'R': 0.22310756972111553}, 967: {'S': 0.5433070866141733}, 975: {'L': 0.5335968379446641}, 980: {'E': 0.5118110236220472}, 983: {'R': 0.24803149606299213}, 985: {'S': 0.3557312252964427}, 986: {'S': 0.24206349206349206}, 989: {'N': 0.37755102040816324}, 992: {'A': 0.3316062176165803}, 993: {'K': 0.39487179487179486}, 994: {'R': 0.3469387755102041}, 995: {'R': 0.40609137055837563}, 1002: {'N': 0.4221105527638191}, 1004: {'D': 0.417910447761194, 'N': 0.373134328358209}, 1006: {'S': 0.42857142857142855}, 1008: {'M': 0.39901477832512317}, 1011: {'R': 0.36764705882352944, 'Q': 0.553921568627451}, 1012: {'A': 0.9509803921568627}, 1018: {'N': 0.3054187192118227}, 1019: {'S': 0.270935960591133}, 1021: {'L': 0.37745098039215685}, 1022: {'G': 0.23383084577114427, 'S': 0.22388059701492538}, 1023: {'A': 0.2885572139303483, 'S': 0.5572139303482587}, 1024: {'E': 0.7892156862745098}, 1025: {'P': 0.4019607843137255}, 1026: {'S': 0.27941176470588236, 'G': 0.25}, 1027: {'F': 0.4068627450980392}, 1028: {'P': 0.33004926108374383}, 1031: {'S': 0.38613861386138615}, 1032: {'S': 0.38308457711442784}, 1034: {'E': 0.43842364532019706}, 1035: {'S': 0.4729064039408867}, 1037: {'V': 0.4079601990049751}, 1040: {'H': 0.26136363636363635}, 1041: {'T': 0.263681592039801, 'Q': 0.3582089552238806}, 1042: {'K': 0.36318407960199006}, 1044: {'M': 0.20833333333333334}, 1046: {'V': 0.35467980295566504}, 1048: {'P': 0.35175879396984927}, 1050: {'A': 0.235, 'M': 0.38}, 1051: {'E': 0.3582089552238806}, 1053: {'S': 0.3448275862068966}, 1055: {'P': 0.35}, 1056: {'S': 0.37438423645320196}, 1057: {'A': 0.3399014778325123}, 1058: {'N': 0.6206896551724138}, 1060: {'R': 0.7303921568627451}, 1061: {'N': 0.36764705882352944}, 1063: {'I': 0.3681592039800995}, 1064: {'G': 0.3910891089108911}, 1065: {'T': 0.22772277227722773}, 1075: {'T': 0.36318407960199006, 'I': 0.25870646766169153}, 1076: {'L': 0.5323383084577115}, 1077: {'H': 0.7810945273631841}}

neutral_clinvar_subs_dict = {250:{'substitution':'K'},
429:{'substitution':'Y'},
445:{'substitution':'A'}, #exist also in msa neutrals
592:{'substitution':'S'},
602:{'substitution':'S'},
942:{'substitution':'K'},
952:{'substitution':'K'},
986:{'substitution':'S'}, #exist also in msa neutrals
996:{'substitution':'S'},
990:{'substitution':'G'},
1000:{'substitution':'G'},
1011:{'substitution':'Q'},#exist also in msa neutrals
1021:{'substitution':'Q'},
1037:{'substitution':'I'},
1027:{'substitution':'I'}}


def create_domain_features():
    topological_domain_dict = {}
    for aa_index in range(1,1079):
        if 1<=aa_index<=558:
            topological_domain_dict[aa_index] = {"VFT":1,"CR":0,"TM":0,"ICL":0,"ECL":0,"Cytoplasmic":0}
        elif 559<=aa_index<=609:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":1,"TM":0,"ICL":0,"ECL":0,"Cytoplasmic":0}
        elif 610<=aa_index<=637:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":1,"ICL":0,"ECL":0,"Cytoplasmic":0}
        elif 638<=aa_index<=646:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":0,"ICL":1,"ECL":0,"Cytoplasmic":0}
        elif 647<=aa_index<=669:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":1,"ICL":0,"ECL":0,"Cytoplasmic":0}
        elif 670<=aa_index<=672:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":0,"ICL":0,"ECL":1,"Cytoplasmic":0}
        elif 673<=aa_index<=706:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":1,"ICL":0,"ECL":0,"Cytoplasmic":0}
        elif 707<=aa_index<=719:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":0,"ICL":1,"ECL":0,"Cytoplasmic":0}
        elif 720<=aa_index<=746:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":1,"ICL":0,"ECL":0,"Cytoplasmic":0}
        elif 747<=aa_index<=769:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":0,"ICL":0,"ECL":1,"Cytoplasmic":0}
        elif 770<=aa_index<=797:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":1,"ICL":0,"ECL":0,"Cytoplasmic":0}
        elif 798<=aa_index<=801:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":0,"ICL":1,"ECL":0,"Cytoplasmic":0}
        elif 802<=aa_index<=828:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":1,"ICL":0,"ECL":0,"Cytoplasmic":0}
        elif 829==aa_index:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":0,"ICL":0,"ECL":1,"Cytoplasmic":0}
        elif 830<=aa_index<=863:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":1,"ICL":0,"ECL":0,"Cytoplasmic":0}
        elif 864<=aa_index<=1078:
            topological_domain_dict[aa_index] = {"VFT":0,"CR":0,"TM":0,"ICL":0,"ECL":0,"Cytoplasmic":1}
    return  topological_domain_dict


blosum62 = substitution_matrices.load(name="BLOSUM62")

def get_fasta_dict(fasta_file):
    fasta_dict = {}
    
    with open(fasta_file,'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.split(">")[1].strip()
                fasta_dict[header] = ''
            else:
                fasta_dict[header] += line.strip()
            
    f.close()

    return fasta_dict



def get_aa_dict(msa_dict,position):
    aa_dict = {}
    for k,v in msa_dict.items():
        aa = v[position]
        if not aa in aa_dict.keys():
            aa_dict[aa] = 1
        else:
            aa_dict[aa] += 1
    if not "-" in aa_dict.keys():
        aa_dict["-"] = 0
    return aa_dict

def compute_identity_score(human_aa_dict,msa_dict,gap_treshold):
    seq_number = len(list(msa_dict.keys()))
    aln_lenght = len(list(msa_dict.values())[0])
    identity_dict = {}
    for i in range(aln_lenght):
        aa_dict = get_aa_dict(msa_dict,i)
        gap_prop = aa_dict["-"]/seq_number
        casr_aa = human_aa_dict[i+1]
        if casr_aa in aa_dict.keys():
            if gap_prop <= gap_treshold:
                identity_score = aa_dict[casr_aa]/(seq_number-aa_dict["-"])
                identity_dict[i+1] = {'identity_score':identity_score,'human_aa':casr_aa}
            elif gap_prop > gap_treshold:
                identity_score = aa_dict[casr_aa]/seq_number
                identity_dict[i+1] = {'identity_score':identity_score,'human_aa':casr_aa}
        else:
            identity_dict[i+1] = {'identity_score':float(0),'human_aa':casr_aa}
    
    return identity_dict


def compute_all_aa_identity_scores(msa_dict,gap_treshold):
    seq_number = len(list(msa_dict.keys()))
    aln_lenght = len(list(msa_dict.values())[0])
    all_aa_identity_dict = {}
    for i in range(aln_lenght):
        aa_dict = get_aa_dict(msa_dict,i)
        gap_prop = aa_dict["-"]/seq_number
        all_aa_identity_dict[i+1] = {}
        for aa, count in aa_dict.items():
            if gap_prop <= gap_treshold:
                prob = count/(seq_number-aa_dict["-"])
                if aa == "-":
                    prob_dict = {aa:float(0)}
                    all_aa_identity_dict[i+1].update(prob_dict)
                else:
                    prob_dict = {aa:prob}
                    all_aa_identity_dict[i+1].update(prob_dict)
            elif gap_prop > gap_treshold:
                prob = count/seq_number
                prob_dict = {aa:prob}
                all_aa_identity_dict[i+1].update(prob_dict)
    return all_aa_identity_dict

def get_substitution_dict(all_aa_identity_dict):
    substitution_dict = {}
    for i in range(1079):
        for pos,consequence in mutations_dict.items():
            position = re.findall(r'\d+', pos)[0]
            substitution = aa_snc_dict[pos[-3:]]
            if i == int(position):
                if substitution in all_aa_identity_dict[i].keys():
                    if i in substitution_dict.keys():
                        substitution_dict[i].update({substitution: all_aa_identity_dict[i][substitution]})
                    else:
                        substitution_dict[i] = {substitution: all_aa_identity_dict[i][substitution]}
                else:
                    if i in substitution_dict.keys():
                        substitution_dict[i].update({substitution:float(0)})
                    else:
                        substitution_dict[i]={substitution:float(0)}

    return substitution_dict
def split_train_validation_positions(tr,fold_seed):
    tr_fold, val_fold = train_test_split(tr, test_size=0.25,random_state=fold_seed,stratify=tr['consequence'])
    return tr_fold, val_fold



def split_train_test_positions(label_pos_csv,seed):
    df = pd.read_csv(label_pos_csv)
    tr, ts = train_test_split(df, test_size=0.2,random_state=seed,stratify=df["consequence"])
    return tr, ts



def split_fasta(fasta_file,seed,fold_seeds):
    fasta_dict = get_fasta_dict(fasta_file)
    key_list = list(fasta_dict.keys())
    df = pd.DataFrame(key_list)
    df_train, df_test = train_test_split(df, test_size=0.2,random_state=seed)

    with open(fasta_file.split(".fasta")[0] + str(seed)+"_test.fasta","w") as f:
        for i in range(len(df_test)) :
            f.write(">" + df_test.iloc[i, 0] + "\n" + fasta_dict[df_test.iloc[i, 0]] + "\n")
    f.close()

    with open(fasta_file.split(".fasta")[0] + str(seed)+"_train.fasta","w") as f:
        for i in range(len(df_train)) :
            f.write(">" + df_train.iloc[i, 0] + "\n" + fasta_dict[df_train.iloc[i, 0]] + "\n")
    f.close()

    for fold_seed in fold_seeds:
        df_train_fold, df_val = train_test_split(df_train, test_size=0.25,random_state=fold_seed)

        with open(fasta_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_train.fasta","w") as f:
            for i in range(len(df_train_fold)) :
                f.write(">" + df_train_fold.iloc[i, 0] + "\n" + fasta_dict[df_train_fold.iloc[i, 0]] + "\n")
        f.close()

        with open(fasta_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_val.fasta","w") as f:
            for i in range(len(df_val)) :
                f.write(">" + df_val.iloc[i, 0] + "\n" + fasta_dict[df_val.iloc[i, 0]] + "\n")
        f.close()    



def write_features_csv(label_pos_csv,casr_file,gprc6a_file,tas1r1_file,tas1r2_file,tas1r3_file,casr_like_file,seed,fold_seed_list):
    topological_domain_dict =  create_domain_features()

    tr, ts = split_train_test_positions(label_pos_csv,seed)

    split_fasta(casr_file,seed,fold_seed_list)
    split_fasta(gprc6a_file,seed,fold_seed_list)
    split_fasta(tas1r1_file,seed,fold_seed_list)
    split_fasta(tas1r2_file,seed,fold_seed_list)
    split_fasta(tas1r3_file,seed,fold_seed_list)
    split_fasta(casr_like_file,seed,fold_seed_list)

    casr_test_dict = get_fasta_dict(casr_file.split(".fasta")[0] + str(seed) +"_test.fasta")
    casr_test_identity_dict = compute_identity_score(human_aa_dict,casr_test_dict,0.3)
    casr_test_all_aa_probs_dict = compute_all_aa_identity_scores(casr_test_dict,0.3)
    casr_test_substiton_dict = get_substitution_dict(casr_test_all_aa_probs_dict)

    gprc6a_test_dict = get_fasta_dict(gprc6a_file.split(".fasta")[0] + str(seed) +"_test.fasta")
    gprc6a_test_identity_dict = compute_identity_score(human_aa_dict,gprc6a_test_dict,0.3)
    gprc6a_test_all_aa_probs_dict = compute_all_aa_identity_scores(gprc6a_test_dict,0.3)
    gprc6a_test_substiton_dict = get_substitution_dict(gprc6a_test_all_aa_probs_dict)

    tas1r1_test_dict = get_fasta_dict(tas1r1_file.split(".fasta")[0] + str(seed) +"_test.fasta")
    tas1r1_test_identity_dict = compute_identity_score(human_aa_dict,tas1r1_test_dict,0.3)
    tas1r1_test_all_aa_probs_dict = compute_all_aa_identity_scores(tas1r1_test_dict,0.3)
    tas1r1_test_substiton_dict = get_substitution_dict(tas1r1_test_all_aa_probs_dict)

    tas1r2_test_dict = get_fasta_dict(tas1r2_file.split(".fasta")[0] + str(seed) +"_test.fasta")
    tas1r2_test_identity_dict = compute_identity_score(human_aa_dict,tas1r2_test_dict,0.3)
    tas1r2_test_all_aa_probs_dict = compute_all_aa_identity_scores(tas1r2_test_dict,0.3)
    tas1r2_test_substiton_dict = get_substitution_dict(tas1r2_test_all_aa_probs_dict)

    tas1r3_test_dict = get_fasta_dict(tas1r3_file.split(".fasta")[0] + str(seed) +"_test.fasta")
    tas1r3_test_identity_dict = compute_identity_score(human_aa_dict,tas1r3_test_dict,0.3)
    tas1r3_test_all_aa_probs_dict = compute_all_aa_identity_scores(tas1r3_test_dict,0.3)
    tas1r3_test_substiton_dict = get_substitution_dict(tas1r3_test_all_aa_probs_dict)

    casr_like_test_dict = get_fasta_dict(casr_like_file.split(".fasta")[0] + str(seed) +"_test.fasta")
    casr_like_test_identity_dict = compute_identity_score(human_aa_dict,casr_like_test_dict,0.3)
    casr_like_test_all_aa_probs_dict = compute_all_aa_identity_scores(casr_like_test_dict,0.3)
    casr_like_test_substiton_dict = get_substitution_dict(casr_like_test_all_aa_probs_dict)

    casr_train_dict = get_fasta_dict(casr_file.split(".fasta")[0] + str(seed) +"_train.fasta")
    casr_train_identity_dict = compute_identity_score(human_aa_dict,casr_train_dict,0.3)
    casr_train_all_aa_probs_dict = compute_all_aa_identity_scores(casr_train_dict,0.3)
    casr_train_substiton_dict = get_substitution_dict(casr_train_all_aa_probs_dict)

    gprc6a_train_dict = get_fasta_dict(gprc6a_file.split(".fasta")[0] + str(seed) +"_train.fasta")
    gprc6a_train_identity_dict = compute_identity_score(human_aa_dict,gprc6a_train_dict,0.3)
    gprc6a_train_all_aa_probs_dict = compute_all_aa_identity_scores(gprc6a_train_dict,0.3)
    gprc6a_train_substiton_dict = get_substitution_dict(gprc6a_train_all_aa_probs_dict)

    tas1r1_train_dict = get_fasta_dict(tas1r1_file.split(".fasta")[0] + str(seed) +"_train.fasta")
    tas1r1_train_identity_dict = compute_identity_score(human_aa_dict,tas1r1_train_dict,0.3)
    tas1r1_train_all_aa_probs_dict = compute_all_aa_identity_scores(tas1r1_train_dict,0.3)
    tas1r1_train_substiton_dict = get_substitution_dict(tas1r1_train_all_aa_probs_dict)

    tas1r2_train_dict = get_fasta_dict(tas1r2_file.split(".fasta")[0] + str(seed) +"_train.fasta")
    tas1r2_train_identity_dict = compute_identity_score(human_aa_dict,tas1r2_train_dict,0.3)
    tas1r2_train_all_aa_probs_dict = compute_all_aa_identity_scores(tas1r2_train_dict,0.3)
    tas1r2_train_substiton_dict = get_substitution_dict(tas1r2_train_all_aa_probs_dict)

    tas1r3_train_dict = get_fasta_dict(tas1r3_file.split(".fasta")[0] + str(seed) +"_train.fasta")
    tas1r3_train_identity_dict = compute_identity_score(human_aa_dict,tas1r3_train_dict,0.3)
    tas1r3_train_all_aa_probs_dict = compute_all_aa_identity_scores(tas1r3_train_dict,0.3)
    tas1r3_train_substiton_dict = get_substitution_dict(tas1r3_train_all_aa_probs_dict)

    casr_like_train_dict = get_fasta_dict(casr_like_file.split(".fasta")[0] + str(seed) +"_train.fasta")
    casr_like_train_identity_dict = compute_identity_score(human_aa_dict,casr_like_train_dict,0.3)
    casr_like_train_all_aa_probs_dict = compute_all_aa_identity_scores(casr_like_train_dict,0.3)
    casr_like_train_substiton_dict = get_substitution_dict(casr_like_train_all_aa_probs_dict)

    for fold_seed in fold_seed_list:
        casr_train_val_dict = get_fasta_dict(casr_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_train.fasta")
        casr_train_val_identity_dict = compute_identity_score(human_aa_dict,casr_train_val_dict,0.3)
        casr_train_val_all_aa_probs_dict = compute_all_aa_identity_scores(casr_train_val_dict,0.3)
        casr_train_val_substiton_dict = get_substitution_dict(casr_train_val_all_aa_probs_dict)

        casr_val_dict = get_fasta_dict(casr_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_val.fasta")
        casr_val_identity_dict = compute_identity_score(human_aa_dict,casr_val_dict,0.3)
        casr_val_all_aa_probs_dict = compute_all_aa_identity_scores(casr_val_dict,0.3)
        casr_val_substiton_dict = get_substitution_dict(casr_val_all_aa_probs_dict)

        gprc6a_train_val_dict = get_fasta_dict(gprc6a_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_train.fasta")
        gprc6a_train_val_identity_dict = compute_identity_score(human_aa_dict,gprc6a_train_val_dict,0.3)
        gprc6a_train_val_all_aa_probs_dict = compute_all_aa_identity_scores(gprc6a_train_val_dict,0.3)
        gprc6a_train_val_substiton_dict = get_substitution_dict(gprc6a_train_val_all_aa_probs_dict)

        gprc6a_val_dict = get_fasta_dict(gprc6a_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_val.fasta")
        gprc6a_val_identity_dict = compute_identity_score(human_aa_dict,gprc6a_val_dict,0.3)
        gprc6a_val_all_aa_probs_dict = compute_all_aa_identity_scores(gprc6a_val_dict,0.3)
        gprc6a_val_substiton_dict = get_substitution_dict(gprc6a_val_all_aa_probs_dict)

        tas1r1_train_val_dict = get_fasta_dict(tas1r1_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_train.fasta")
        tas1r1_train_val_identity_dict = compute_identity_score(human_aa_dict,tas1r1_train_val_dict,0.3)
        tas1r1_train_val_all_aa_probs_dict = compute_all_aa_identity_scores(tas1r1_train_val_dict,0.3)
        tas1r1_train_val_substiton_dict = get_substitution_dict(tas1r1_train_val_all_aa_probs_dict)

        tas1r1_val_dict = get_fasta_dict(tas1r1_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_val.fasta")
        tas1r1_val_identity_dict = compute_identity_score(human_aa_dict,tas1r1_val_dict,0.3)
        tas1r1_val_all_aa_probs_dict = compute_all_aa_identity_scores(tas1r1_val_dict,0.3)
        tas1r1_val_substiton_dict = get_substitution_dict(tas1r1_val_all_aa_probs_dict)

      

        tas1r2_train_val_dict = get_fasta_dict(tas1r2_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_train.fasta")
        tas1r2_train_val_identity_dict = compute_identity_score(human_aa_dict,tas1r2_train_val_dict,0.3)
        tas1r2_train_val_all_aa_probs_dict = compute_all_aa_identity_scores(tas1r2_train_val_dict,0.3)
        tas1r2_train_val_substiton_dict = get_substitution_dict(tas1r2_train_val_all_aa_probs_dict)

        tas1r2_val_dict = get_fasta_dict(tas1r2_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_val.fasta")
        tas1r2_val_identity_dict = compute_identity_score(human_aa_dict,tas1r2_val_dict,0.3)
        tas1r2_val_all_aa_probs_dict = compute_all_aa_identity_scores(tas1r2_val_dict,0.3)
        tas1r2_val_substiton_dict = get_substitution_dict(tas1r2_val_all_aa_probs_dict)

        

        tas1r3_train_val_dict = get_fasta_dict(tas1r3_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_train.fasta")
        tas1r3_train_val_identity_dict = compute_identity_score(human_aa_dict,tas1r3_train_val_dict,0.3)
        tas1r3_train_val_all_aa_probs_dict = compute_all_aa_identity_scores(tas1r3_train_val_dict,0.3)
        tas1r3_train_val_substiton_dict = get_substitution_dict(tas1r3_train_val_all_aa_probs_dict)

        tas1r3_val_dict = get_fasta_dict(tas1r3_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_val.fasta")
        tas1r3_val_identity_dict = compute_identity_score(human_aa_dict,tas1r3_val_dict,0.3)
        tas1r3_val_all_aa_probs_dict = compute_all_aa_identity_scores(tas1r3_val_dict,0.3)
        tas1r3_val_substiton_dict = get_substitution_dict(tas1r3_val_all_aa_probs_dict)

        

        casr_like_train_val_dict = get_fasta_dict(casr_like_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_train.fasta")
        casr_like_train_val_identity_dict = compute_identity_score(human_aa_dict,casr_like_train_val_dict,0.3)
        casr_like_train_val_all_aa_probs_dict = compute_all_aa_identity_scores(casr_like_train_val_dict,0.3)
        casr_like_train_val_substiton_dict = get_substitution_dict(casr_like_train_val_all_aa_probs_dict)

        casr_like_val_dict = get_fasta_dict(casr_like_file.split(".fasta")[0] + str(seed)+"_"+str(fold_seed) +"_val.fasta")
        casr_like_val_identity_dict = compute_identity_score(human_aa_dict,casr_like_val_dict,0.3)
        casr_like_val_all_aa_probs_dict = compute_all_aa_identity_scores(casr_like_val_dict,0.3)
        casr_like_val_substiton_dict = get_substitution_dict(casr_like_val_all_aa_probs_dict)

        tr_fold, val_fold = split_train_validation_positions(tr,fold_seed)
        columns = ["position","substitution","consequence"]
        tr_fold_dict_list = tr_fold[columns].to_dict(orient='records')
        tr_val_fold_list = val_fold[columns].to_dict(orient='records')

        with open("/cta/users/abircan/casr_artcile_ml_check_11_06_2023/train_csv_files/"+str(seed)+"_"+str(fold_seed)+"_train.csv",'w') as f:
            casr_pred_writer = csv.writer(f, delimiter=',')
            casr_pred_writer = csv.writer(f, delimiter=',')
            casr_pred_writer.writerow(["position","aa_0_A",	"aa_0_R","aa_0_N", "aa_0_D", "aa_0_C", "aa_0_Q","aa_0_E", "aa_0_G",	"aa_0_H","aa_0_I","aa_0_L",	"aa_0_K","aa_0_M","aa_0_F",
            "aa_0_P","aa_0_S","aa_0_T","aa_0_W","aa_0_Y","aa_0_V","aa_0_B","aa_0_Z","aa_0_X","aa_0_*",
            "CaSR_a0","GPRC6A_a0","TAS1R1_a0","TAS1R2_a0","TAS1R3_a0","CaSR_like_a0",
            "zimmerman_polarity_a0","average_flexibility_a0",
            "dayhoff_a0", "avg_buried_area_a0","hydropathicity_a0","atomic_weight_ratio_a0",
            "molecular_weight_a0","bulkiness_a0","tm_tendency_a0",
            "VFT", "CR","TM","ICL","ECL","Cytoplasmic",
            "aa_1_A",	"aa_1_R","aa_1_N", "aa_1_D", "aa_1_C", "aa_1_Q","aa_1_E", "aa_1_G",	"aa_1_H","aa_1_I","aa_1_L",	"aa_1_K","aa_1_M","aa_1_F",
            "aa_1_P","aa_1_S","aa_1_T","aa_1_W","aa_1_Y","aa_1_V","aa_1_B","aa_1_Z","aa_1_X","aa_1_*",
            "CaSR_a1","GPRC6A_a1","TAS1R1_a1","TAS1R2_a1","TAS1R3_a1","CaSR_like_a1",
            "zimmerman_polarity_a1","average_flexibility_a1",
            "dayhoff_a1", "avg_buried_area_a1","hydropathicity_a1","atomic_weight_ratio_a1",
            "molecular_weight_a1","bulkiness_a1","tm_tendency_a1",
            "consequence"])
            for _ in tr_fold_dict_list:
                for pos,consequence in mutations_dict.items():  
                    position = re.findall(r'\d+', pos)[0]
                    substitution = aa_snc_dict[pos[-3:]]
                    if int(position) == _['position'] and substitution ==  _['substitution'] and disease_cause_dict[consequence] ==  _['consequence']:
                        casr_pred_writer.writerow([int(position),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[0]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[1]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[2]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[3]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[4]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[5]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[6]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[7]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[8]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[9]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[10]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[11]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[12]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[13]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[14]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[15]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[16]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[17]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[18]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[19]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[20]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[21]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[22]),
                        float(blosum62.__getitem__(casr_train_val_identity_dict[int(position)]['human_aa']).values()[23]),



                        casr_train_val_identity_dict[int(position)]['identity_score'],
                        gprc6a_train_val_identity_dict[int(position)]['identity_score'], tas1r1_train_val_identity_dict[int(position)]['identity_score'],
                        tas1r2_train_val_identity_dict[int(position)]['identity_score'],tas1r3_train_val_identity_dict[int(position)]['identity_score'],
                        casr_like_train_val_identity_dict[int(position)]['identity_score'],
                        zimmerman_polarity_dict[casr_train_val_identity_dict[int(position)]['human_aa']],
                        average_flexibility_dict[casr_train_val_identity_dict[int(position)]['human_aa']],
                        dayhoff_dict[casr_train_val_identity_dict[int(position)]['human_aa']],
                        avg_buried_area_dict[casr_train_val_identity_dict[int(position)]['human_aa']],
                        doolittle_hydropathicity_dict[casr_train_val_identity_dict[int(position)]['human_aa']],
                        atomic_weight_ratio_dict[casr_train_val_identity_dict[int(position)]['human_aa']],
                        molecular_weight_dict[casr_train_val_identity_dict[int(position)]['human_aa']],
                        bulkiness_dict[casr_train_val_identity_dict[int(position)]['human_aa']],
                        tm_tendency_dict[casr_train_val_identity_dict[int(position)]['human_aa']],

                        topological_domain_dict[int(position)]["VFT"], topological_domain_dict[int(position)]["CR"],
                        topological_domain_dict[int(position)]["TM"], topological_domain_dict[int(position)]["ICL"],
                        topological_domain_dict[int(position)]["ECL"], topological_domain_dict[int(position)]["Cytoplasmic"],

                        float(blosum62.__getitem__(substitution).values()[0]),
                        float(blosum62.__getitem__(substitution).values()[1]),
                        float(blosum62.__getitem__(substitution).values()[2]),
                        float(blosum62.__getitem__(substitution).values()[3]),
                        float(blosum62.__getitem__(substitution).values()[4]),
                        float(blosum62.__getitem__(substitution).values()[5]),
                        float(blosum62.__getitem__(substitution).values()[6]),
                        float(blosum62.__getitem__(substitution).values()[7]),
                        float(blosum62.__getitem__(substitution).values()[8]),
                        float(blosum62.__getitem__(substitution).values()[9]),
                        float(blosum62.__getitem__(substitution).values()[10]),
                        float(blosum62.__getitem__(substitution).values()[11]),
                        float(blosum62.__getitem__(substitution).values()[12]),
                        float(blosum62.__getitem__(substitution).values()[13]),
                        float(blosum62.__getitem__(substitution).values()[14]),
                        float(blosum62.__getitem__(substitution).values()[15]),
                        float(blosum62.__getitem__(substitution).values()[16]),
                        float(blosum62.__getitem__(substitution).values()[17]),
                        float(blosum62.__getitem__(substitution).values()[18]),
                        float(blosum62.__getitem__(substitution).values()[19]),
                        float(blosum62.__getitem__(substitution).values()[20]),
                        float(blosum62.__getitem__(substitution).values()[21]),
                        float(blosum62.__getitem__(substitution).values()[22]),
                        float(blosum62.__getitem__(substitution).values()[23]),
                        
                        casr_train_val_substiton_dict[int(position)][substitution],
                        gprc6a_train_val_substiton_dict[int(position)][substitution],
                        tas1r1_train_val_substiton_dict[int(position)][substitution],
                        tas1r2_train_val_substiton_dict[int(position)][substitution],
                        tas1r3_train_val_substiton_dict[int(position)][substitution],
                        casr_like_train_val_substiton_dict[int(position)][substitution],

                        zimmerman_polarity_dict[substitution],
                        average_flexibility_dict[substitution],
                        dayhoff_dict[substitution],
                        avg_buried_area_dict[substitution],
                        doolittle_hydropathicity_dict[substitution],
                        atomic_weight_ratio_dict[substitution],
                        molecular_weight_dict[substitution],
                        bulkiness_dict[substitution],
                        tm_tendency_dict[substitution],
                        disease_cause_dict[consequence]]) 

        f.close()
        with open("/cta/users/abircan/casr_artcile_ml_check_11_06_2023/train_csv_files/"+str(seed)+"_"+str(fold_seed)+"_val.csv",'w') as f:
            casr_pred_writer = csv.writer(f, delimiter=',')
            casr_pred_writer.writerow(["position","aa_0_A",	"aa_0_R","aa_0_N", "aa_0_D", "aa_0_C", "aa_0_Q","aa_0_E", "aa_0_G",	"aa_0_H","aa_0_I","aa_0_L",	"aa_0_K","aa_0_M","aa_0_F",
            "aa_0_P","aa_0_S","aa_0_T","aa_0_W","aa_0_Y","aa_0_V","aa_0_B","aa_0_Z","aa_0_X","aa_0_*",
            "CaSR_a0","GPRC6A_a0","TAS1R1_a0","TAS1R2_a0","TAS1R3_a0","CaSR_like_a0",
            "zimmerman_polarity_a0","average_flexibility_a0",
            "dayhoff_a0", "avg_buried_area_a0","hydropathicity_a0","atomic_weight_ratio_a0",
            "molecular_weight_a0","bulkiness_a0","tm_tendency_a0",
            "VFT", "CR","TM","ICL","ECL","Cytoplasmic",
            "aa_1_A",	"aa_1_R","aa_1_N", "aa_1_D", "aa_1_C", "aa_1_Q","aa_1_E", "aa_1_G",	"aa_1_H","aa_1_I","aa_1_L",	"aa_1_K","aa_1_M","aa_1_F",
            "aa_1_P","aa_1_S","aa_1_T","aa_1_W","aa_1_Y","aa_1_V","aa_1_B","aa_1_Z","aa_1_X","aa_1_*",
            "CaSR_a1","GPRC6A_a1","TAS1R1_a1","TAS1R2_a1","TAS1R3_a1","CaSR_like_a1",
            "zimmerman_polarity_a1","average_flexibility_a1",
            "dayhoff_a1", "avg_buried_area_a1","hydropathicity_a1","atomic_weight_ratio_a1",
            "molecular_weight_a1","bulkiness_a1","tm_tendency_a1",
            "consequence"])
            for _ in tr_val_fold_list:
                for pos,consequence in mutations_dict.items():  
                    position = re.findall(r'\d+', pos)[0]
                    substitution = aa_snc_dict[pos[-3:]]
                    if int(position) == _['position'] and substitution ==  _['substitution'] and disease_cause_dict[consequence] ==  _['consequence']:
                        casr_pred_writer.writerow([int(position),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[0]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[1]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[2]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[3]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[4]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[5]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[6]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[7]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[8]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[9]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[10]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[11]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[12]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[13]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[14]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[15]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[16]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[17]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[18]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[19]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[20]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[21]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[22]),
                        float(blosum62.__getitem__(casr_val_identity_dict[int(position)]['human_aa']).values()[23]),



                        casr_val_identity_dict[int(position)]['identity_score'],
                        gprc6a_val_identity_dict[int(position)]['identity_score'], tas1r1_val_identity_dict[int(position)]['identity_score'],
                        tas1r2_val_identity_dict[int(position)]['identity_score'],tas1r3_val_identity_dict[int(position)]['identity_score'],
                        casr_like_val_identity_dict[int(position)]['identity_score'],
                        zimmerman_polarity_dict[casr_val_identity_dict[int(position)]['human_aa']],
                        average_flexibility_dict[casr_val_identity_dict[int(position)]['human_aa']],
                        dayhoff_dict[casr_val_identity_dict[int(position)]['human_aa']],
                        avg_buried_area_dict[casr_val_identity_dict[int(position)]['human_aa']],
                        doolittle_hydropathicity_dict[casr_val_identity_dict[int(position)]['human_aa']],
                        atomic_weight_ratio_dict[casr_val_identity_dict[int(position)]['human_aa']],
                        molecular_weight_dict[casr_val_identity_dict[int(position)]['human_aa']],
                        bulkiness_dict[casr_val_identity_dict[int(position)]['human_aa']],
                        tm_tendency_dict[casr_val_identity_dict[int(position)]['human_aa']],
                        topological_domain_dict[int(position)]["VFT"], topological_domain_dict[int(position)]["CR"],
                        topological_domain_dict[int(position)]["TM"], topological_domain_dict[int(position)]["ICL"],
                        topological_domain_dict[int(position)]["ECL"], topological_domain_dict[int(position)]["Cytoplasmic"],
    
                        float(blosum62.__getitem__(substitution).values()[0]),
                        float(blosum62.__getitem__(substitution).values()[1]),
                        float(blosum62.__getitem__(substitution).values()[2]),
                        float(blosum62.__getitem__(substitution).values()[3]),
                        float(blosum62.__getitem__(substitution).values()[4]),
                        float(blosum62.__getitem__(substitution).values()[5]),
                        float(blosum62.__getitem__(substitution).values()[6]),
                        float(blosum62.__getitem__(substitution).values()[7]),
                        float(blosum62.__getitem__(substitution).values()[8]),
                        float(blosum62.__getitem__(substitution).values()[9]),
                        float(blosum62.__getitem__(substitution).values()[10]),
                        float(blosum62.__getitem__(substitution).values()[11]),
                        float(blosum62.__getitem__(substitution).values()[12]),
                        float(blosum62.__getitem__(substitution).values()[13]),
                        float(blosum62.__getitem__(substitution).values()[14]),
                        float(blosum62.__getitem__(substitution).values()[15]),
                        float(blosum62.__getitem__(substitution).values()[16]),
                        float(blosum62.__getitem__(substitution).values()[17]),
                        float(blosum62.__getitem__(substitution).values()[18]),
                        float(blosum62.__getitem__(substitution).values()[19]),
                        float(blosum62.__getitem__(substitution).values()[20]),
                        float(blosum62.__getitem__(substitution).values()[21]),
                        float(blosum62.__getitem__(substitution).values()[22]),
                        float(blosum62.__getitem__(substitution).values()[23]),
                        
                        casr_val_substiton_dict[int(position)][substitution],
                        gprc6a_val_substiton_dict[int(position)][substitution],
                        tas1r1_val_substiton_dict[int(position)][substitution],
                        tas1r2_val_substiton_dict[int(position)][substitution],
                        tas1r3_val_substiton_dict[int(position)][substitution],
                        casr_like_val_substiton_dict[int(position)][substitution],

                        zimmerman_polarity_dict[substitution],
                        average_flexibility_dict[substitution],
                        dayhoff_dict[substitution],
                        avg_buried_area_dict[substitution],
                        doolittle_hydropathicity_dict[substitution],
                        atomic_weight_ratio_dict[substitution],
                        molecular_weight_dict[substitution],
                        bulkiness_dict[substitution],
                        tm_tendency_dict[substitution],
                        disease_cause_dict[consequence]]) 

        f.close()
    tr_dict_list = tr[columns].to_dict(orient='records')
    ts_dict_list = ts[columns].to_dict(orient='records')
    with open("/cta/users/abircan/casr_artcile_ml_check_11_06_2023/train_csv_files/"+str(seed)+"_test.csv",'w') as f:
        casr_pred_writer = csv.writer(f, delimiter=',')
        casr_pred_writer.writerow(["position","aa_0_A",	"aa_0_R","aa_0_N", "aa_0_D", "aa_0_C", "aa_0_Q","aa_0_E", "aa_0_G",	"aa_0_H","aa_0_I","aa_0_L",	"aa_0_K","aa_0_M","aa_0_F",
        "aa_0_P","aa_0_S","aa_0_T","aa_0_W","aa_0_Y","aa_0_V","aa_0_B","aa_0_Z","aa_0_X","aa_0_*",
        "CaSR_a0","GPRC6A_a0","TAS1R1_a0","TAS1R2_a0","TAS1R3_a0","CaSR_like_a0",
        "zimmerman_polarity_a0","average_flexibility_a0",
        "dayhoff_a0", "avg_buried_area_a0","hydropathicity_a0","atomic_weight_ratio_a0",
        "molecular_weight_a0","bulkiness_a0","tm_tendency_a0",
        "VFT", "CR","TM","ICL","ECL","Cytoplasmic",
        "aa_1_A",	"aa_1_R","aa_1_N", "aa_1_D", "aa_1_C", "aa_1_Q","aa_1_E", "aa_1_G",	"aa_1_H","aa_1_I","aa_1_L",	"aa_1_K","aa_1_M","aa_1_F",
        "aa_1_P","aa_1_S","aa_1_T","aa_1_W","aa_1_Y","aa_1_V","aa_1_B","aa_1_Z","aa_1_X","aa_1_*",
        "CaSR_a1","GPRC6A_a1","TAS1R1_a1","TAS1R2_a1","TAS1R3_a1","CaSR_like_a1",
        "zimmerman_polarity_a1","average_flexibility_a1",
        "dayhoff_a1", "avg_buried_area_a1","hydropathicity_a1","atomic_weight_ratio_a1",
        "molecular_weight_a1","bulkiness_a1","tm_tendency_a1",
        "consequence"])
        for _ in ts_dict_list:
            for pos,consequence in mutations_dict.items():  
                position = re.findall(r'\d+', pos)[0]
                substitution = aa_snc_dict[pos[-3:]]
                if int(position) == _['position'] and substitution ==  _['substitution'] and disease_cause_dict[consequence] ==  _['consequence']:
                    casr_pred_writer.writerow([int(position),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[0]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[1]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[2]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[3]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[4]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[5]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[6]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[7]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[8]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[9]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[10]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[11]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[12]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[13]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[14]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[15]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[16]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[17]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[18]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[19]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[20]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[21]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[22]),
                    float(blosum62.__getitem__(casr_test_identity_dict[int(position)]['human_aa']).values()[23]),
                
                    casr_test_identity_dict[int(position)]['identity_score'],
                    gprc6a_test_identity_dict[int(position)]['identity_score'], tas1r1_test_identity_dict[int(position)]['identity_score'],
                    tas1r2_test_identity_dict[int(position)]['identity_score'],tas1r3_test_identity_dict[int(position)]['identity_score'],
                    casr_like_test_identity_dict[int(position)]['identity_score'],
                    zimmerman_polarity_dict[casr_test_identity_dict[int(position)]['human_aa']],
                    average_flexibility_dict[casr_test_identity_dict[int(position)]['human_aa']],
                    dayhoff_dict[casr_test_identity_dict[int(position)]['human_aa']],
                    avg_buried_area_dict[casr_test_identity_dict[int(position)]['human_aa']],
                    doolittle_hydropathicity_dict[casr_test_identity_dict[int(position)]['human_aa']],
                    atomic_weight_ratio_dict[casr_test_identity_dict[int(position)]['human_aa']],
                    molecular_weight_dict[casr_test_identity_dict[int(position)]['human_aa']],
                    bulkiness_dict[casr_test_identity_dict[int(position)]['human_aa']],
                    tm_tendency_dict[casr_test_identity_dict[int(position)]['human_aa']],
                    
                    topological_domain_dict[int(position)]["VFT"], topological_domain_dict[int(position)]["CR"],
                    topological_domain_dict[int(position)]["TM"], topological_domain_dict[int(position)]["ICL"],
                    topological_domain_dict[int(position)]["ECL"], topological_domain_dict[int(position)]["Cytoplasmic"],

                    float(blosum62.__getitem__(substitution).values()[0]),
                    float(blosum62.__getitem__(substitution).values()[1]),
                    float(blosum62.__getitem__(substitution).values()[2]),
                    float(blosum62.__getitem__(substitution).values()[3]),
                    float(blosum62.__getitem__(substitution).values()[4]),
                    float(blosum62.__getitem__(substitution).values()[5]),
                    float(blosum62.__getitem__(substitution).values()[6]),
                    float(blosum62.__getitem__(substitution).values()[7]),
                    float(blosum62.__getitem__(substitution).values()[8]),
                    float(blosum62.__getitem__(substitution).values()[9]),
                    float(blosum62.__getitem__(substitution).values()[10]),
                    float(blosum62.__getitem__(substitution).values()[11]),
                    float(blosum62.__getitem__(substitution).values()[12]),
                    float(blosum62.__getitem__(substitution).values()[13]),
                    float(blosum62.__getitem__(substitution).values()[14]),
                    float(blosum62.__getitem__(substitution).values()[15]),
                    float(blosum62.__getitem__(substitution).values()[16]),
                    float(blosum62.__getitem__(substitution).values()[17]),
                    float(blosum62.__getitem__(substitution).values()[18]),
                    float(blosum62.__getitem__(substitution).values()[19]),
                    float(blosum62.__getitem__(substitution).values()[20]),
                    float(blosum62.__getitem__(substitution).values()[21]),
                    float(blosum62.__getitem__(substitution).values()[22]),
                    float(blosum62.__getitem__(substitution).values()[23]),
                    
                    casr_test_substiton_dict[int(position)][substitution],
                    gprc6a_test_substiton_dict[int(position)][substitution],
                    tas1r1_test_substiton_dict[int(position)][substitution],
                    tas1r2_test_substiton_dict[int(position)][substitution],
                    tas1r3_test_substiton_dict[int(position)][substitution],
                    casr_like_test_substiton_dict[int(position)][substitution],

                    zimmerman_polarity_dict[substitution],
                    average_flexibility_dict[substitution],
                    dayhoff_dict[substitution],
                    avg_buried_area_dict[substitution],
                    doolittle_hydropathicity_dict[substitution],
                    atomic_weight_ratio_dict[substitution],
                    molecular_weight_dict[substitution],
                    bulkiness_dict[substitution],
                    tm_tendency_dict[substitution],
                    disease_cause_dict[consequence]]) 

    f.close()
    with open("/cta/users/abircan/casr_artcile_ml_check_11_06_2023/train_csv_files/"+str(seed)+"_train.csv",'w') as f:
        casr_pred_writer = csv.writer(f, delimiter=',')
        casr_pred_writer.writerow(["position","aa_0_A",	"aa_0_R","aa_0_N", "aa_0_D", "aa_0_C", "aa_0_Q","aa_0_E", "aa_0_G",	"aa_0_H","aa_0_I","aa_0_L",	"aa_0_K","aa_0_M","aa_0_F",
        "aa_0_P","aa_0_S","aa_0_T","aa_0_W","aa_0_Y","aa_0_V","aa_0_B","aa_0_Z","aa_0_X","aa_0_*",
        "CaSR_a0","GPRC6A_a0","TAS1R1_a0","TAS1R2_a0","TAS1R3_a0","CaSR_like_a0",
        "zimmerman_polarity_a0","average_flexibility_a0",
        "dayhoff_a0", "avg_buried_area_a0","hydropathicity_a0","atomic_weight_ratio_a0",
        "molecular_weight_a0","bulkiness_a0","tm_tendency_a0",
        "VFT", "CR","TM","ICL","ECL","Cytoplasmic",
        "aa_1_A",	"aa_1_R","aa_1_N", "aa_1_D", "aa_1_C", "aa_1_Q","aa_1_E", "aa_1_G",	"aa_1_H","aa_1_I","aa_1_L",	"aa_1_K","aa_1_M","aa_1_F",
        "aa_1_P","aa_1_S","aa_1_T","aa_1_W","aa_1_Y","aa_1_V","aa_1_B","aa_1_Z","aa_1_X","aa_1_*",
        "CaSR_a1","GPRC6A_a1","TAS1R1_a1","TAS1R2_a1","TAS1R3_a1","CaSR_like_a1",
        "zimmerman_polarity_a1","average_flexibility_a1",
        "dayhoff_a1", "avg_buried_area_a1","hydropathicity_a1","atomic_weight_ratio_a1",
        "molecular_weight_a1","bulkiness_a1","tm_tendency_a1",
        "consequence"])
        for _ in tr_dict_list:  
            for pos,consequence in mutations_dict.items():  
                position = re.findall(r'\d+', pos)[0]
                substitution = aa_snc_dict[pos[-3:]]
                if int(position) == _['position'] and substitution ==  _['substitution'] and disease_cause_dict[consequence] ==  _['consequence']:
                    casr_pred_writer.writerow([int(position),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[0]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[1]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[2]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[3]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[4]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[5]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[6]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[7]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[8]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[9]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[10]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[11]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[12]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[13]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[14]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[15]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[16]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[17]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[18]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[19]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[20]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[21]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[22]),
                    float(blosum62.__getitem__(casr_train_identity_dict[int(position)]['human_aa']).values()[23]),



                    casr_train_identity_dict[int(position)]['identity_score'],
                    gprc6a_train_identity_dict[int(position)]['identity_score'], tas1r1_train_identity_dict[int(position)]['identity_score'],
                    tas1r2_train_identity_dict[int(position)]['identity_score'],tas1r3_train_identity_dict[int(position)]['identity_score'],
                    casr_like_train_identity_dict[int(position)]['identity_score'],
                    zimmerman_polarity_dict[casr_train_identity_dict[int(position)]['human_aa']],
                    average_flexibility_dict[casr_train_identity_dict[int(position)]['human_aa']],
                    dayhoff_dict[casr_train_identity_dict[int(position)]['human_aa']],
                    avg_buried_area_dict[casr_train_identity_dict[int(position)]['human_aa']],
                    doolittle_hydropathicity_dict[casr_train_identity_dict[int(position)]['human_aa']],
                    atomic_weight_ratio_dict[casr_train_identity_dict[int(position)]['human_aa']],
                    molecular_weight_dict[casr_train_identity_dict[int(position)]['human_aa']],
                    bulkiness_dict[casr_train_identity_dict[int(position)]['human_aa']],
                    tm_tendency_dict[casr_train_identity_dict[int(position)]['human_aa']],

                    topological_domain_dict[int(position)]["VFT"], topological_domain_dict[int(position)]["CR"],
                    topological_domain_dict[int(position)]["TM"], topological_domain_dict[int(position)]["ICL"],
                    topological_domain_dict[int(position)]["ECL"], topological_domain_dict[int(position)]["Cytoplasmic"],
 
                    float(blosum62.__getitem__(substitution).values()[0]),
                    float(blosum62.__getitem__(substitution).values()[1]),
                    float(blosum62.__getitem__(substitution).values()[2]),
                    float(blosum62.__getitem__(substitution).values()[3]),
                    float(blosum62.__getitem__(substitution).values()[4]),
                    float(blosum62.__getitem__(substitution).values()[5]),
                    float(blosum62.__getitem__(substitution).values()[6]),
                    float(blosum62.__getitem__(substitution).values()[7]),
                    float(blosum62.__getitem__(substitution).values()[8]),
                    float(blosum62.__getitem__(substitution).values()[9]),
                    float(blosum62.__getitem__(substitution).values()[10]),
                    float(blosum62.__getitem__(substitution).values()[11]),
                    float(blosum62.__getitem__(substitution).values()[12]),
                    float(blosum62.__getitem__(substitution).values()[13]),
                    float(blosum62.__getitem__(substitution).values()[14]),
                    float(blosum62.__getitem__(substitution).values()[15]),
                    float(blosum62.__getitem__(substitution).values()[16]),
                    float(blosum62.__getitem__(substitution).values()[17]),
                    float(blosum62.__getitem__(substitution).values()[18]),
                    float(blosum62.__getitem__(substitution).values()[19]),
                    float(blosum62.__getitem__(substitution).values()[20]),
                    float(blosum62.__getitem__(substitution).values()[21]),
                    float(blosum62.__getitem__(substitution).values()[22]),
                    float(blosum62.__getitem__(substitution).values()[23]),
                    
                    casr_train_substiton_dict[int(position)][substitution],
                    gprc6a_train_substiton_dict[int(position)][substitution],
                    tas1r1_train_substiton_dict[int(position)][substitution],
                    tas1r2_train_substiton_dict[int(position)][substitution],
                    tas1r3_train_substiton_dict[int(position)][substitution],
                    casr_like_train_substiton_dict[int(position)][substitution],

                    zimmerman_polarity_dict[substitution],
                    average_flexibility_dict[substitution],
                    dayhoff_dict[substitution],
                    avg_buried_area_dict[substitution],
                    doolittle_hydropathicity_dict[substitution],
                    atomic_weight_ratio_dict[substitution],
                    molecular_weight_dict[substitution],
                    bulkiness_dict[substitution],
                    tm_tendency_dict[substitution],
                    disease_cause_dict[consequence]]) 

    f.close()




def generate_random_seed(reperat_number):
    seed_list = []
    for i in range(1,reperat_number+1):
        seed_list.append(1000+i)
    return seed_list

seed_list = generate_random_seed(50)

for seed in seed_list:
    write_features_csv("/cta/users/abircan/casr_artcile_ml_check_11_06_2023/position_label.csv","/cta/users/abircan/casr_artcile_ml_check_11_06_2023/CaSR_nogapcasr.fasta","/cta/users/abircan/casr_artcile_ml_check_11_06_2023/GPRC6A_nogapcasr.fasta","/cta/users/abircan/casr_artcile_ml_check_11_06_2023/TAS1R1_nogapcasr.fasta","/cta/users/abircan/casr_artcile_ml_check_11_06_2023/TAS1R2_nogapcasr.fasta","/cta/users/abircan/casr_artcile_ml_check_11_06_2023/TAS1R3_nogapcasr.fasta","/cta/users/abircan/casr_artcile_ml_check_11_06_2023/all_casr_likes_nogapcasr.fasta",seed,[1,2,3,4,5])



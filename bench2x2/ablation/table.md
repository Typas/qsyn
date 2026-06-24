# Gadgetize merge-strategy ablation (T-count primary; lowest wins)

| circuit    | budget  | strategy      |  t_count |    peak_kb |     wall | status |
|------------|---------|---------------|----------|------------|----------|--------|
| hwb8       | 256M    | z-overlap     |     3420 |     169080 |  1:33.94 | ok |
| hwb8       | 256M    | max-merge     |        ? |     207664 | 30:00.01 | TLE |
| hwb8       | 256M    | max-terms     |        ? |     228084 | 30:00.01 | TLE |
| hwb8       | 256M    | balanced-fill |     3505 |     162800 |  0:48.04 | ok |
| ham15-high | 64M     | z-overlap     |      949 |      37584 |  4:04.05 | ok |
| ham15-high | 64M     | max-merge     |      823 |      58408 | 16:47.52 | ok |
| ham15-high | 64M     | max-terms     |        ? |      70488 | 30:00.00 | TLE |
| ham15-high | 64M     | balanced-fill |     1010 |      36212 |  0:08.23 | ok |

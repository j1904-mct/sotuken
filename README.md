# sotuken

基本的にはJun_Makino/treecode.Cの203行目にあるcalculate_gravity()の実行時間を計測する。
ただし、1回目のmallocなどに時間を掛けないように、何回かcalculate_gravity()を実行する。

calculate_gravity()周辺で時間計測を行うようにして(<chrono>ヘッダを利用)、Nを同様に変えたときにキャッシュを利用しているときと利用していないときの時間差を考える。
その結果からキャッシュブロッキングした時の性能の予測モデルを作成できるとよい。

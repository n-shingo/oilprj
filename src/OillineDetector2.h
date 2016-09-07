#ifndef __OILLINE_DETECTOR2_H__
#define __OILLINE_DETECTOR2_H__

#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

namespace sn{


	// エッジ線検出器クラス
	class OillineDetector2{

		/////////////////
		// 公開メソッド  //
		/////////////////
	public:

		// コンストラクタ
		OillineDetector2();

		// エッジ抽出実行
		int Execute(Mat &src, double *dist, double *gl_theta, Mat &result_img, bool debug=false);

		/////////////////
		// プロパティ   //
		/////////////////
		
		// 検出する最大スリットの数
		void SetMaxSlitCount( int cnt ){ _max_slit_cnt= cnt; }
		int GetMaxSlitCount( void ){ return _max_slit_cnt; }


		// 俯瞰画像のプロパティ
		void SetCamTopForBird(int top){ _camTop = top; set_coeffs(); }  // 俯瞰前の画像の俯瞰画像上部y座標
		void SetBirdHeight(int height){ _birdH = height; set_coeffs(); }  // 俯瞰画像の高さ[pix]
		void SetBirdBtmX(int left, int right){ _birdBtmLft = left; _birdBtmRgt = right; set_coeffs(); }  // 俯瞰画像下部の左右x座標

		// 距離変換のプロパティ
		void SetDgl(double d_gl){ _d_gl = d_gl; }  // D_gl(車軸とカメラ画像最下部までの距離[mm])を設定
		void SetDppX(double dpp_x){ _dpp_x = dpp_x; } // DppX(x軸方向の距離変換係数[mm/pix])を設定
		void SetDppY(double dpp_y){ _dpp_y = dpp_y; } // DppY(y軸方向の距離変換係数[mm/pix])を設定


	private:
		////////////////////
		//  非公開メソッド  //
		////////////////////

		// ２つの画像が同じ内容かチェックする
		bool is_same_image( Mat &img1, Mat &img2 );

		// 必要な係数などを設定する
		void set_coeffs();
		
		// ガンマ補正を行う
		Mat auto_gamma_correction(Mat &src, double base, double scale);
		
		// 位置xの縦ラインピークをcnt個取得する．ただしth以上の値.
		vector<int> find_vertical_peaks( Mat &gray, int x, int minInterval, int cnt, int th );
		
		// 距離th以下のピーク値をチェーン化する
		// margin, step は入力するピーク値の最左座標[pix]と間隔[pix]
		vector< vector<cv::Point> > make_chains( vector< vector<int> > &peaks, int th, int margin, int step );

		// ぐちゃぐちゃしたチェーンをchainsから削除する．[長さ/端点直線距離]がthより大きいと削除される.
		void remove_noisy_chains( vector< vector<Point> > &chains, double th );
		
		// 二値化処理によってのっぺりとしたピークを削除する（線でないものを削除する）
		// 閾値は二値化した時の黒ピクセルの割合で、thより大きいとチェーンからピークが削除される
		void remove_flat_peaks( vector< vector<Point> > &chains, Mat &gray, double th );
		
		// 平均値をしきい値とする二値化
		void binarize_by_average( Mat &gray, int x, int y, int size, Mat &dst );

		// 平均輝度の高いcnt個のチェーンに絞る
		void focus_top_chains( vector< vector<cv::Point> > &chains, Mat &gray, int cnt );

		// 線形補間により間を補完する
		vector< vector<double> > interpolate_chains( vector< vector<Point> > &chains );

		// 平滑化をする
		vector<double> smoothing( vector<double> &data, int win_size );

		// インターバル間での曲率を求める
		vector<double> calc_curvature(vector<double> ps, int interval);

		// エッジ線が許容範囲内で動いたかチェック
		bool is_acceptable_difference( Vec2d l1, Vec2d l2, double th1, double th2 );

		// カメラ画像の歪み補正をする(0:失敗、1:成功)
		void cameara_calibrate(Mat& src, Mat& dst);

		// 俯瞰画像を作成する(0:失敗、1:成功)
		void make_birdimg(Mat& src, Mat& dst);

		// 元座標から俯瞰画像へ座標変換する
		vector<Point2d> to_bird_coordinate(vector<Point2d> &points);
		
		// 画像上の直線のrho(>0)とtheta(0-2PI)を求める
		bool get_rho_theta( vector<Point2d> &pnts, double* rho, double* theta );

		// 画像上の線から、実空間上での線の値と距離を求める
		// rho, theta :画像上の線,  rho_gl, theta_gl: 実空間上の線(結果), dist:距離(結果)
		// 戻り値 0:該当なし、1:成功
		void real_world_value(double rho, double theta, double* rho_gl, double* theta_gl, double* dist);

		// 直線を描画する
		Mat draw_lines(Mat img, vector<Vec2d> lines, Scalar color, int thickness = 1);




		////////////////////
		//   メンバ変数    //
		////////////////////
		
		// 最大検出スリット数
		int _max_slit_cnt;

		// 最後に処理した画像
		Mat _last_img;

		// 最後の実行結果
		bool _last_exec_ret;
		double _last_dist;
		double _last_gl_theta;
		Mat _last_result_img;
			
		// 最後に検出された線
		Vec2d _last_line;

		// ガンマ補正
		double _gamma_base; // ガンマ補正基準値
		double _gamma_scale; // ガンマ補正スケール

		// ガウシアン
		int _gaussian_window;  // 窓サイズ(奇数)
		double _gaussian_sigma; // 標準偏差

		// 縦ラインのピーク検出
		int _peak_margin;   // ピーク検出の最左点座標
		int _peak_interval; // ピーク検出間隔
		int _peak_isolate_interval; // ピーク間の最低距離
		int _peak_count;    // 検出するピークの数
		int _peak_min;      // ピークの最低輝度値 

		// チェーン化
		int _chain_max_interval; // チェーン化の最高距離
		int _chain_min_size;     // チェーンの最低サイズ
		double _chain_noisy_th;  // ぐちゃぐちゃしたチェーンの閾値
		
		// ピークの削除
		double _flat_peak_th;    // のっぺりピークの閾値

		// 二値化
		int _binarize_size;  // 局所的二値化を適用する窓サイズ[pix]

		// 平滑化
		int _smooth_win_size; // 窓サイズ
		
		// 曲率
		int _curve_interval; // 曲率のためのインターバル
		double _th_curve_max; // 溶接溝における曲率の最低値

		// 前回線との変化許容閾値
		double _th_rho_mm;
		double _th_theta_deg;


		// カメラ校正パラメータ
		Mat _calibIntrinsic;
		Mat _calibDistortion;

		// カメラ画像サイズ
		int _camW;
		int _camH;

		// 俯瞰変換
		int _camTop;     // 俯瞰前の上部y座標
		int _birdW;      // 俯瞰画像幅
		int _birdH;      // 俯瞰画像高さ
		int _birdBtmLft; // 俯瞰画像左下座標
		int _birdBtmRgt; // 俯瞰画像右下座標
		Mat _homoMat;    // ホモグラフィ行列

		// 距離算出に関する変数
		double _d_gl;  // 車軸とカメラ画像最下部までの距離[mm]
		double _dpp_x; // x軸方向の距離変換係数[mm/pix]
		double _dpp_y; // y軸方向の距離変換係数[mm/pix]
	};

} // nemespace gp


#endif // __OILLINE_DTECTOR2_H__

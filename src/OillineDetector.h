#ifndef __OILLINE_DETECTOR_H__
#define __OILLINE_DETECTOR_H__

#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

class OillineDetector{

	/////////////////
	// 公開メソッド  //
	/////////////////
public:

	// コンストラクタ
	OillineDetector();

	// エッジ抽出実行
	int Execute(Mat &src, double *dist, double *gl_theta, Mat &result_img);

	/////////////////
	// プロパティ   //
	/////////////////

	// 俯瞰画像のプロパティ
	void SetCamTopForBird(int top){ _camTop = top; set_coeffs(); }  // 俯瞰前の画像の俯瞰画像上部y座標
	void SetBirdHeight(int height){ _birdH = height; set_coeffs(); }  // 俯瞰画像の高さ[pix]
	void SetBirdBtmX(int left, int right){ _birdBtmLft = left; _birdBtmRgt = right; set_coeffs(); }  // 俯瞰画像下部の左右x座標

	// 距離変換のプロパティ
	void SetDgl(double d_gl){ _d_gl = d_gl; }  // D_gl(車軸とカメラ画像最下部までの距離[mm])を設定
	void SetDppX(double dpp_x){ _dpp_x = dpp_x; } // DppX(x軸方向の距離変換係数[mm/pix])を設定
	void SetDppY(double dpp_y){ _dpp_y = dpp_y; } // DppY(y軸方向の距離変換係数[mm/pix])を設定

	// 画像下部の無視する領域高さ[pix]
	void SetIgnoreBottomSize( int size ){ _ignore_bottom_size = size; }
	int GetIgnoreBottomSize( void ){ return _ignore_bottom_size; }
	
	// ガンマ補正 プロパティ
	void SetGammaBase( double base ){ _gamma_base = base; }
	double GetGammaBase( void ){ return _gamma_base; }
	void SetGammaScale( double scale ){ _gamma_scale = scale; }
	double GetGammaScale( void ){ return _gamma_scale; }
	
	// Gaussianフィルタ プロパティ
	void SetGaussianWindow(int win){ _gaussian_window = win; }
	int GetGaussianWindow(void){ return _gaussian_window; }
	void SetGaussianSigma(double sigma){ _gaussian_sigma = sigma; }
	double GetGaussianSigma(void){ return _gaussian_sigma; }


private:
	////////////////////
	//  非公開メソッド  //
	////////////////////

	// 必要な係数などを設定する
	void set_coeffs();

	// カメラ画像の歪み補正をする(0:失敗、1:成功)
	void cameara_calibrate(Mat& src, Mat& dst);

	// 俯瞰画像を作成する(0:失敗、1:成功)
	void make_birdimg(Mat& src, Mat& dst);

	// 元座標から俯瞰画像へ座標変換する
	vector<Point2d> to_bird_coordinate(vector<Point2d> &points);
	
	// 画像上の直線のrho(>0)とtheta(0-2PI)を求める
	void get_rho_theta(Point2d p1, Point2d p2, double *rho, double *theta);

	// 画像上の線から、実空間上での線の値と距離を求める
	// rho, theta :画像上の線,  rho_gl, theta_gl: 実空間上の線(結果), dist:距離(結果)
	// 戻り値 0:該当なし、1:成功
	void real_world_value(double rho, double theta, double* rho_gl, double* theta_gl, double* dist);

	// ガンマ補正を行う
	Mat auto_gamma_correction(Mat &src, double base, double scale);

	// hsvの閾値処理によってスリットラインを抽出した8bit画像を返す
	// Hue[0-179], Sat[0-255], Val[0-255]において, th_hue_low < Hue <= th_hue_low && Sat <= th_sat && Val <= th_val
	Mat hsv_slitline_threshold(Mat &src, int th_hue_low, int th_hue_up, int th_sat, int th_val);
	Mat hsv_slitline_threshold2(Mat &src, int th_hue_low, int th_hue_up, int th_sat, int th_val);

	// つながっている領域毎に画像を分割する.
	// dstは分割された画像群, rectsは元画像に対する矩形領域を表す
	// ドットのような1点のデータは出力データに含めない
	void split_image_region(Mat &gray, vector< Mat > &dst, vector< Rect > &rects);

	// 領域内の黒領域を塗りつぶす
	void fill_region(Mat &gray);

	// 細線化（水平方向に移動しつつ平均位置をプロットする）
	vector<double> thining_horizontal_line(Mat &gray);

	// インターバル間での曲率を求める
	vector<double> calc_curvature(vector<double> ps, int interval);

	// 曲率のグラフを画像に描画する
	Mat draw_curvature_plot(Mat img, vector<double> curvature, Scalar line_color, double threshold = CV_PI, Scalar th_color = Scalar(0, 0, 255));

	// 分割された画像を合成する. 黒色部分は合成の対象から外す．
	Mat fusion_splited_images(Size image_size, vector<Mat> images, vector<Rect> rects);

	// 直線を描画する
	Mat draw_lines(Mat img, vector<Vec2f> lines, Scalar color, int thickness = 1);




	////////////////////
	//   メンバ変数    //
	////////////////////

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

	// 処理を無視する下部の画像サイズ
	int _ignore_bottom_size;

	// ガンマ補正
	double _gamma_base; // ガンマ補正基準値
	double _gamma_scale; // ガンマ補正スケール

	// ガウシアン
	int _gaussian_window;  // 窓サイズ(奇数)
	double _gaussian_sigma; // 標準偏差

	// HSV閾値
	int _th_hue_low;  // hue下限
	int _th_hue_up;   // hue上限
	int _th_sat;      // sat下限
	int _th_val;      // val下限

	// 分割した領域の最低サイズ
	int _min_region_size;

	// 曲率
	int _curve_interval; // 曲率のためのインターバル
	double _th_curve_max; // 溶接溝における曲率の最低値


	// 距離算出に関する変数
	double _d_gl;  // 車軸とカメラ画像最下部までの距離[mm]
	double _dpp_x; // x軸方向の距離変換係数[mm/pix]
	double _dpp_y; // y軸方向の距離変換係数[mm/pix]


};



#endif // __OILLINE_DTECTOR_H__

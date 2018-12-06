/*
    Jon Sweeney
    November 8, 2018
    DM : EM implementation for both discrete and continuous data sets.
         -Does not have a general solution for all data sets but it sufficient for the assignment.
*/

#include <iostream> 
#include <vector> 
#include <cmath>

using namespace std;

// Contains the relative frequency of heads vs tails for each set of 10 coin flips.
struct CoinFlip {
    float prob_heads;
    float prob_tails;
};

struct CoinFlipResults {
    float num_heads;
    float num_tails;
};

struct MaximizationTuple {
    CoinFlipResults coinA;
    CoinFlipResults coinB;
};

struct MeanTuple {
    float mean_one;
    float mean_two;
};

// Need to iterate till difference between terms is < 0.01, switching this from 0.1 to allow for more iterations.
const float EPSILON = 0.01;
const vector<float> weather_data = {70, 62, 89, 54, 97, 75, 82, 56, 32,78};
const vector<CoinFlip> test_data = { {0.5, 0.5}, {0.9, 0.1}, {0.8, 0.2}, {0.4, 0.6}, {0.7, 0.3}};
//-------------------Begin functionality for problem 1---------------------------------------

inline float binomialDist(float prob, int n, int x) {
    return pow(prob, x) * pow((1 - prob), (n - x));
}

/* Need to do this for all 5 sets of 10 coin flips each iteration. */
void probExpectation(float A_PROB, float B_PROB, int n, int x, vector<MaximizationTuple> &max_set) {
    // First calculate the next value for A_PROB
    float new_prob_A = binomialDist(A_PROB, n, x);

    // Now calculate the next value for B_PROB
    float new_prob_B = binomialDist(B_PROB, n, x);

    // Normalize the new probabilities.
    float normalization = new_prob_A + new_prob_B;
    new_prob_A /= normalization;
    new_prob_B /= normalization;

    // Now estimate the *likely* # heads vs # tails for coin A and coin B
    float a_heads = new_prob_A * x;
    float a_tails = new_prob_A * (n - x);
    CoinFlipResults a_flip = {a_heads, a_tails};

    // Further build the new set of expected heads and tails for each coin.
    float b_heads = new_prob_B * x;
    float b_tails = new_prob_B * (n - x);
    CoinFlipResults b_flip = {b_heads, b_tails};

    //NOTE: This builds the set of num_heads vs num_tails for coins A and B,
    // During maximization I need to tally these up and get the new A_PROB and B_PROB
    max_set.push_back({a_flip, b_flip});
}

/* Now get the new A_PROB and B_PROB based off the maximization_set */
void probMaximization(float &A_PROB, float &B_PROB, vector<MaximizationTuple> &max_set) {
    float num_aHeads = 0;
    float num_aTails = 0;
    float num_bHeads = 0;
    float num_bTails = 0;

    for (auto itr: max_set) {
        num_aHeads += itr.coinA.num_heads;
        num_aTails += itr.coinA.num_tails;
        num_bHeads += itr.coinB.num_heads;
        num_bTails += itr.coinB.num_tails;
    }

    // Now set up the new A_PROB and B_PROB
    A_PROB = num_aHeads / (num_aHeads + num_aTails);
    B_PROB = num_bHeads / (num_bHeads + num_bTails);

    max_set.clear();
}

//-------------------Begin functionality for problem 2---------------------------------------

inline float numerator(float mean, float std_dev, float x) {
    return (exp((-1.0 / (2.0 * (pow(std_dev, 2.0)))) * pow(x - mean, 2.0)));
}

void expectation(vector<float> &init_means, vector<float> &sun_ex, vector<float> &clo_ex, float std_dev, int k) {
    float num = 0.0;
    float denom = 0.0;
    // First compute expectation set for a sunny day.
    for (auto itr: weather_data) {
        num = 0.0;
        denom = 0.0;
        num = numerator(init_means[0], std_dev, itr);
        for (int i = 0; i < k; i++) {
            // Now get the denominator, not a good general solution, but makes it easier to achieve.
            denom += (exp((-1.0 / (2.0 * (pow(std_dev, 2.0)))) * pow(itr - init_means[i], 2.0)));
        }
        float test = num / denom;

        sun_ex.push_back(test);
    }

    // Now compute it for a cloudy day.

    for (auto itr: weather_data) {
        num = 0.0;
        denom = 0.0;
        num = numerator(init_means[1], std_dev, itr);
        //Now get denom
        for (int i = 0; i < k; i++) {
            denom += (exp((-1.0 / (2.0 * (pow(std_dev, 2.0)))) * pow(itr - init_means[i], 2.0)));
        }
        float test = num / denom;

        clo_ex.push_back(test);
    }
}

void maximization(vector<float> &init_means, vector<float> &sun_ex, vector<float> &clo_ex, vector<MeanTuple> &mean_list) {
    // Iterate through the sunny expectation set and assign a new mean.
    float sun_total = 0.0;
    float sun_mean = 0.0;
    for (int i = 0; i < sun_ex.size(); i++) {
        sun_total += sun_ex[i];
        sun_mean += (sun_ex[i] * weather_data[i]);
    }
    sun_mean /= sun_total;

    // Now do the same for cloudy expectation set.
    float clo_total = 0.0;
    float clo_mean = 0.0;

    for (int i = 0; i < clo_ex.size(); i++) {
        clo_total += clo_ex[i];
        clo_mean += (clo_ex[i] * weather_data[i]);
    }
    clo_mean /= clo_total;

    mean_list.push_back({sun_mean, clo_mean});

}

int main() {

    // Begin work for first problem
    vector<MaximizationTuple> max_set = {};
    // The initial assignments for probability of heads on coins A and B
    float A_PROB = 0.6;
    float B_PROB = 0.5;
    vector<CoinFlipResults> coin_prob;

    for (int i = 0; i < 10; i++) {
        for (auto itr: test_data) {
            probExpectation(A_PROB, B_PROB, 10, itr.prob_heads * 10, max_set);
        }
        probMaximization(A_PROB, B_PROB, max_set);
        coin_prob.push_back({A_PROB, B_PROB});
    }
    cout << "Coin probability of heads." << endl;
    cout << "A_PROB(Heads)" << "\t\t" << "B_PROB(Heads)" << endl;
    for (auto itr: coin_prob) {
        cout << itr.num_heads << "\t\t" << itr.num_tails << endl;
    }

    // Begin work for second problem.
    vector<float> sunny_expectation = {};
    vector<float> cloudy_expectation = {};
    vector<MeanTuple> mean_list = {};
    vector<float> init_means = {80.0, 55.0};

    // Set up my initial guesses for the means.
    mean_list.push_back({80.0, 55.0});

    float sun_control = 0.0;
    float clo_control = 0.0;

    while ((fabs(mean_list.back().mean_one - sun_control) > EPSILON) && (fabs(mean_list.back().mean_two - clo_control) > EPSILON)) {
        // Set up my loop control.
        sun_control = mean_list.back().mean_one;
        clo_control = mean_list.back().mean_two;
        // Ensure my expectation sets are empty after previous iterations
        sunny_expectation.clear();
        cloudy_expectation.clear();

        // First build expectation sets
        expectation(init_means, sunny_expectation, cloudy_expectation, 10, 2);

        // Now get the new means
        maximization(init_means, sunny_expectation, cloudy_expectation, mean_list);

        init_means.clear();
        init_means.push_back(mean_list.back().mean_one);
        init_means.push_back(mean_list.back().mean_two);
    }

    // Check the resulting mean of each iteration now.
    cout << "Sunny " << "\t\t" << "Cloudy" << endl;
    for (auto itr: mean_list) {
        cout << itr.mean_one << "\t\t" << itr.mean_two << endl;
    }

}

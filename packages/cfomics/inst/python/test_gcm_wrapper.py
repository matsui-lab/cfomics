import pytest
import numpy as np
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gcm_dowhy import (
    run_dowhy_gcm,
    validate_inputs,
    load_dag,
    load_config,
    process_X,
    summarize_ite_distribution,
    bootstrap_ate
)


@pytest.fixture
def sample_data():
    """Generate sample data for testing."""
    np.random.seed(42)
    n = 100
    X = np.random.randn(n, 3)
    T = (X[:, 0] > 0).astype(int)
    Y = 1.5 * T + 0.5 * X[:, 1] + np.random.randn(n) * 0.1
    return X, T, Y


@pytest.fixture
def sample_dag_file(tmp_path):
    """Create a temporary DAG file."""
    dag_path = tmp_path / "test_dag.json"
    dag_config = {
        "nodes": ["X0", "X1", "X2", "T", "Y"],
        "edges": [["X0", "T"], ["T", "Y"], ["X1", "Y"]]
    }
    with open(dag_path, 'w') as f:
        json.dump(dag_config, f)
    return str(dag_path)


@pytest.fixture
def sample_config_file(tmp_path):
    """Create a temporary config file."""
    config_path = tmp_path / "test_config.json"
    config = {
        "random_seed": 42,
        "bootstrap": {
            "n_iterations": 10
        },
        "dimensionality_reduction": {
            "enabled": False
        },
        "output_metadata": False
    }
    with open(config_path, 'w') as f:
        json.dump(config, f)
    return str(config_path)


class TestValidateInputs:
    """Tests for input validation function."""
    
    def test_valid_inputs(self, sample_data):
        """Test validation with valid inputs."""
        X, T, Y = sample_data
        X_out, T_out, Y_out = validate_inputs(X, T, Y)
        assert X_out.shape == X.shape
        assert T_out.shape == T.shape
        assert Y_out.shape == Y.shape
    
    def test_shape_mismatch_T(self, sample_data):
        """Test validation catches T shape mismatch."""
        X, T, Y = sample_data
        with pytest.raises(ValueError, match="Shape mismatch"):
            validate_inputs(X, T[:-1], Y)
    
    def test_shape_mismatch_Y(self, sample_data):
        """Test validation catches Y shape mismatch."""
        X, T, Y = sample_data
        with pytest.raises(ValueError, match="Shape mismatch"):
            validate_inputs(X, T, Y[:-1])
    
    def test_missing_values_X(self, sample_data):
        """Test validation catches missing values in X."""
        X, T, Y = sample_data
        X_nan = X.copy()
        X_nan[0, 0] = np.nan
        with pytest.raises(ValueError, match="missing values"):
            validate_inputs(X_nan, T, Y)
    
    def test_missing_values_T(self, sample_data):
        """Test validation catches missing values in T."""
        X, T, Y = sample_data
        T_nan = T.astype(float)
        T_nan[0] = np.nan
        with pytest.raises(ValueError, match="missing values"):
            validate_inputs(X, T_nan, Y)
    
    def test_missing_values_Y(self, sample_data):
        """Test validation catches missing values in Y."""
        X, T, Y = sample_data
        Y_nan = Y.copy()
        Y_nan[0] = np.nan
        with pytest.raises(ValueError, match="missing values"):
            validate_inputs(X, T, Y_nan)
    
    def test_non_binary_treatment(self, sample_data):
        """Test validation catches non-binary treatment."""
        X, T, Y = sample_data
        T_bad = T.copy().astype(float)
        T_bad[0] = 2
        with pytest.raises(ValueError, match="binary"):
            validate_inputs(X, T_bad, Y)
    
    def test_converts_lists_to_arrays(self):
        """Test that lists are converted to numpy arrays."""
        X = [[1, 2], [3, 4], [5, 6]]
        T = [0, 1, 0]
        Y = [1.0, 2.0, 1.5]
        X_out, T_out, Y_out = validate_inputs(X, T, Y)
        assert isinstance(X_out, np.ndarray)
        assert isinstance(T_out, np.ndarray)
        assert isinstance(Y_out, np.ndarray)
    
    def test_1d_X_reshaped(self):
        """Test that 1D X is reshaped to 2D."""
        X = np.array([1, 2, 3])
        T = np.array([0, 1, 0])
        Y = np.array([1.0, 2.0, 1.5])
        X_out, _, _ = validate_inputs(X, T, Y)
        assert X_out.ndim == 2
        assert X_out.shape == (3, 1)


class TestLoadDag:
    """Tests for DAG loading function."""
    
    def test_load_dag_success(self, sample_dag_file):
        """Test successful DAG loading from JSON."""
        graph = load_dag(sample_dag_file)
        assert len(graph.nodes()) == 5
        assert len(graph.edges()) == 3
        assert ("X0", "T") in graph.edges()
        assert ("T", "Y") in graph.edges()
        assert ("X1", "Y") in graph.edges()
    
    def test_load_dag_file_not_found(self):
        """Test error when DAG file doesn't exist."""
        with pytest.raises(FileNotFoundError):
            load_dag("/nonexistent/path/dag.json")
    
    def test_load_dag_invalid_format(self, tmp_path):
        """Test error when DAG file has invalid format."""
        invalid_dag = tmp_path / "invalid_dag.json"
        with open(invalid_dag, 'w') as f:
            json.dump({"nodes": ["A", "B"]}, f)
        with pytest.raises(ValueError, match="edges"):
            load_dag(str(invalid_dag))


class TestLoadConfig:
    """Tests for config loading function."""
    
    def test_load_config_default(self):
        """Test loading default config when no file provided."""
        config = load_config(None)
        assert 'random_seed' in config
        assert 'bootstrap' in config
        assert 'dimensionality_reduction' in config
    
    def test_load_config_from_file(self, sample_config_file):
        """Test loading config from file."""
        config = load_config(sample_config_file)
        assert config['random_seed'] == 42
        assert config['bootstrap']['n_iterations'] == 10
    
    def test_load_config_nonexistent_file(self):
        """Test that nonexistent file returns default config."""
        config = load_config("/nonexistent/config.json")
        assert 'random_seed' in config


class TestProcessX:
    """Tests for X processing function."""
    
    def test_no_reduction(self):
        """Test X processing without dimensionality reduction."""
        X = np.random.randn(100, 5)
        config = {'dimensionality_reduction': {'enabled': False}}
        X_processed, metadata = process_X(X, config)
        
        np.testing.assert_array_equal(X, X_processed)
        assert metadata['original_shape'] == (100, 5)
    
    def test_pca_reduction(self):
        """Test X processing with PCA dimensionality reduction."""
        X = np.random.randn(100, 20)
        config = {
            'dimensionality_reduction': {
                'enabled': True,
                'method': 'pca',
                'n_components': 5
            }
        }
        X_processed, metadata = process_X(X, config)
        
        assert X_processed.shape == (100, 5)
        assert metadata['method'] == 'pca'
        assert 'explained_variance_ratio' in metadata
        assert len(metadata['explained_variance_ratio']) == 5
    
    def test_pca_components_capped(self):
        """Test that PCA components are capped at min(n_samples, n_features)."""
        X = np.random.randn(50, 10)
        config = {
            'dimensionality_reduction': {
                'enabled': True,
                'method': 'pca',
                'n_components': 100
            }
        }
        X_processed, metadata = process_X(X, config)
        
        assert X_processed.shape[1] <= min(50, 10)
    
    def test_autoencoder_not_implemented(self):
        """Test that autoencoder raises NotImplementedError."""
        X = np.random.randn(100, 20)
        config = {
            'dimensionality_reduction': {
                'enabled': True,
                'method': 'autoencoder',
                'n_components': 5
            }
        }
        with pytest.raises(NotImplementedError):
            process_X(X, config)
    
    def test_unknown_method(self):
        """Test that unknown method raises ValueError."""
        X = np.random.randn(100, 20)
        config = {
            'dimensionality_reduction': {
                'enabled': True,
                'method': 'unknown_method',
                'n_components': 5
            }
        }
        with pytest.raises(ValueError, match="Unknown"):
            process_X(X, config)


class TestSummarizeIteDistribution:
    """Tests for ITE distribution summary function."""
    
    def test_summary_keys(self):
        """Test that summary contains all expected keys."""
        ite = np.random.randn(1000)
        summary = summarize_ite_distribution(ite)
        
        expected_keys = ['mean', 'std', 'quantile_05', 'quantile_50', 'quantile_95', 'min', 'max']
        for key in expected_keys:
            assert key in summary
    
    def test_summary_values_reasonable(self):
        """Test that summary values are reasonable for normal distribution."""
        np.random.seed(42)
        ite = np.random.randn(10000)
        summary = summarize_ite_distribution(ite)
        
        assert abs(summary['mean']) < 0.1
        assert 0.9 < summary['std'] < 1.1
        assert abs(summary['quantile_50'] - summary['mean']) < 0.1
    
    def test_summary_handles_list(self):
        """Test that summary handles list input."""
        ite = [1.0, 2.0, 3.0, 4.0, 5.0]
        summary = summarize_ite_distribution(ite)
        assert summary['mean'] == 3.0


class TestRunDowhyGcm:
    """Tests for main run_dowhy_gcm function."""
    
    def test_basic_execution(self, sample_data):
        """Test basic execution of run_dowhy_gcm."""
        X, T, Y = sample_data
        result = run_dowhy_gcm(X, T, Y, return_metadata=False)
        
        assert "ite" in result
        assert "ate" in result
        assert "y0_hat" in result
        assert "y1_hat" in result
        assert "summary" in result
        assert "samples" in result
        
        assert len(result["ite"]) == len(Y)
        assert isinstance(result["ate"], float)
    
    def test_summary_structure(self, sample_data):
        """Test that summary has correct structure."""
        X, T, Y = sample_data
        result = run_dowhy_gcm(X, T, Y, return_metadata=False)
        
        summary = result["summary"]
        assert "ate" in summary
        assert "ate_ci_lower" in summary
        assert "ate_ci_upper" in summary
        assert "ite_mean" in summary
        assert "ite_std" in summary
        assert "ite_quantiles" in summary
    
    def test_samples_structure(self, sample_data):
        """Test that samples has correct structure."""
        X, T, Y = sample_data
        result = run_dowhy_gcm(X, T, Y, return_metadata=False)
        
        samples = result["samples"]
        assert "y0" in samples
        assert "y1" in samples
        assert "ite" in samples
        assert len(samples["y0"]) == len(Y)
    
    def test_metadata_included(self, sample_data):
        """Test that metadata is included when requested."""
        X, T, Y = sample_data
        result = run_dowhy_gcm(X, T, Y, return_metadata=True)
        
        assert "metadata" in result
        metadata = result["metadata"]
        assert "input_hash" in metadata
        assert "dag_structure" in metadata
        assert "random_seed" in metadata
    
    def test_deterministic_with_seed(self, sample_data, sample_config_file):
        """Test that results are deterministic with same seed."""
        X, T, Y = sample_data
        
        result1 = run_dowhy_gcm(X, T, Y, config_path=sample_config_file, return_metadata=False)
        result2 = run_dowhy_gcm(X, T, Y, config_path=sample_config_file, return_metadata=False)
        
        np.testing.assert_array_almost_equal(
            result1["ite"],
            result2["ite"],
            decimal=5
        )
    
    def test_ate_reasonable_range(self, sample_data):
        """Test that ATE is in reasonable range for synthetic data."""
        X, T, Y = sample_data
        result = run_dowhy_gcm(X, T, Y, return_metadata=False)
        
        ate = result["summary"]["ate"]
        assert 0.5 < ate < 2.5, f"ATE {ate} outside expected range [0.5, 2.5]"
    
    def test_with_custom_dag(self, sample_data, sample_dag_file):
        """Test execution with custom DAG file."""
        X, T, Y = sample_data
        result = run_dowhy_gcm(X, T, Y, dag_path=sample_dag_file, return_metadata=True)
        
        assert "metadata" in result
        assert "dag_structure" in result["metadata"]
    
    def test_backward_compatibility(self, sample_data):
        """Test backward compatibility with original signature."""
        X, T, Y = sample_data
        result = run_dowhy_gcm(X, T, Y)
        
        assert "ite" in result
        assert "ate" in result
        assert "y0_hat" in result
        assert "y1_hat" in result
        
        assert isinstance(result["ite"], np.ndarray)
        assert isinstance(result["ate"], float)


class TestVisualization:
    """Tests for visualization functions."""
    
    def test_plot_ite_distribution(self, sample_data, tmp_path):
        """Test ITE distribution plot."""
        from visualization import plot_ite_distribution
        
        X, T, Y = sample_data
        ite = np.random.randn(100)
        
        save_path = str(tmp_path / "ite_dist.png")
        fig, ax = plot_ite_distribution(ite, save_path=save_path)
        
        assert fig is not None
        assert ax is not None
        assert os.path.exists(save_path)
        
        plt_module = __import__('matplotlib.pyplot', fromlist=[''])
        plt_module.close(fig)
    
    def test_plot_outcome_shift(self, tmp_path):
        """Test outcome shift plot."""
        from visualization import plot_outcome_shift
        
        y0 = np.random.randn(100)
        y1 = y0 + 1.5 + np.random.randn(100) * 0.5
        
        save_path = str(tmp_path / "outcome_shift.png")
        fig, axes = plot_outcome_shift(y0, y1, save_path=save_path)
        
        assert fig is not None
        assert len(axes) == 2
        assert os.path.exists(save_path)
        
        plt_module = __import__('matplotlib.pyplot', fromlist=[''])
        plt_module.close(fig)
    
    def test_plot_uplift_curve(self, tmp_path):
        """Test uplift curve plot."""
        from visualization import plot_uplift_curve
        
        ite = np.random.randn(100)
        
        save_path = str(tmp_path / "uplift.png")
        fig, ax = plot_uplift_curve(ite, save_path=save_path)
        
        assert fig is not None
        assert ax is not None
        assert os.path.exists(save_path)
        
        plt_module = __import__('matplotlib.pyplot', fromlist=[''])
        plt_module.close(fig)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

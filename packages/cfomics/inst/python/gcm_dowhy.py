# inst/python/gcm_dowhy.py
# Enhanced DoWhy GCM wrapper with bootstrap, PCA, and metadata support
import numpy as np
import pandas as pd
import networkx as nx
from dowhy import gcm
import json
import os
import hashlib
from datetime import datetime
import inspect


def _supports_rng_parameter():
    """Check if gcm.interventional_samples supports the rng parameter."""
    try:
        sig = inspect.signature(gcm.interventional_samples)
        return 'rng' in sig.parameters
    except Exception:
        return False


def _interventional_samples_compat(scm, interventions, num_samples_to_draw, random_state=None):
    """
    Version-compatible wrapper for gcm.interventional_samples.
    
    Handles both newer DoWhy versions (with rng parameter) and older versions
    (without rng parameter, using global numpy seed).
    """
    if _supports_rng_parameter() and random_state is not None:
        rng = np.random.default_rng(random_state)
        return gcm.interventional_samples(
            scm, interventions, num_samples_to_draw=num_samples_to_draw, rng=rng
        )
    else:
        if random_state is not None:
            np.random.seed(random_state)
        return gcm.interventional_samples(
            scm, interventions, num_samples_to_draw=num_samples_to_draw
        )


def load_config(path=None):
    """
    Load configuration from JSON file.
    
    Parameters
    ----------
    path : str, optional
        Path to JSON config file. If None, returns default config.
        
    Returns
    -------
    dict
        Configuration dictionary
    """
    default_config = {
        'dag_file': None,
        'interventions': {
            'type': 'binary',
            'values': [0, 1],
            'treatment_variable': 'T'
        },
        'bootstrap': {
            'enabled': False,
            'n_iterations': 100,
            'confidence_level': 0.95
        },
        'dimensionality_reduction': {
            'enabled': False,
            'method': 'pca',
            'n_components': 10
        },
        'random_seed': 0,
        'output_metadata': False,
        'metadata_output_path': 'metadata.json'
    }
    
    if path is None or not os.path.exists(path):
        return default_config
    
    with open(path, 'r') as f:
        user_config = json.load(f)
    
    for key, value in user_config.items():
        if isinstance(value, dict) and key in default_config:
            default_config[key].update(value)
        else:
            default_config[key] = value
    
    return default_config


def load_dag(path):
    """
    Load DAG structure from JSON file.
    
    Parameters
    ----------
    path : str
        Path to JSON file containing DAG specification.
        Expected format: {"nodes": [...], "edges": [[from, to], ...]}
        
    Returns
    -------
    networkx.DiGraph
        The causal graph
        
    Raises
    ------
    FileNotFoundError
        If the DAG file does not exist
    ValueError
        If the DAG file format is invalid
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"DAG file not found: {path}")
    
    with open(path, 'r') as f:
        dag_config = json.load(f)
    
    if 'edges' not in dag_config:
        raise ValueError("DAG config must contain 'edges' key")
    
    graph = nx.DiGraph()
    
    if 'nodes' in dag_config:
        graph.add_nodes_from(dag_config['nodes'])
    
    graph.add_edges_from(dag_config['edges'])
    
    return graph


def validate_inputs(X, T, Y):
    """
    Validate input data for causal inference.
    
    Parameters
    ----------
    X : numpy.ndarray
        Covariates matrix (n, p)
    T : numpy.ndarray
        Treatment vector (n,)
    Y : numpy.ndarray
        Outcome vector (n,)
        
    Returns
    -------
    tuple
        Validated (X, T, Y) as numpy arrays
        
    Raises
    ------
    ValueError
        If inputs have invalid shapes, contain missing values, or T is not binary
    """
    if not isinstance(X, np.ndarray):
        X = np.asarray(X)
    if not isinstance(T, np.ndarray):
        T = np.asarray(T)
    if not isinstance(Y, np.ndarray):
        Y = np.asarray(Y)
    
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    
    if X.ndim != 2:
        raise ValueError(f"X must be 2D, got shape {X.shape}")
    if T.ndim != 1:
        raise ValueError(f"T must be 1D, got shape {T.shape}")
    if Y.ndim != 1:
        raise ValueError(f"Y must be 1D, got shape {Y.shape}")
    
    n = X.shape[0]
    if T.shape[0] != n:
        raise ValueError(f"Shape mismatch: X has {n} samples but T has {T.shape[0]}")
    if Y.shape[0] != n:
        raise ValueError(f"Shape mismatch: X has {n} samples but Y has {Y.shape[0]}")
    
    if np.isnan(X).any():
        raise ValueError("X contains missing values (NaN)")
    if np.isnan(T).any():
        raise ValueError("T contains missing values (NaN)")
    if np.isnan(Y).any():
        raise ValueError("Y contains missing values (NaN)")
    
    unique_T = np.unique(T)
    if not np.all(np.isin(unique_T, [0, 1])):
        raise ValueError(f"T must be binary (0/1), got unique values: {unique_T}")
    
    return np.asarray(X), np.asarray(T), np.asarray(Y)


def process_X(X, config):
    """
    Process X matrix with optional dimensionality reduction.
    
    Parameters
    ----------
    X : numpy.ndarray
        Input covariate matrix (n, p)
    config : dict
        Configuration dictionary with 'dimensionality_reduction' settings
        
    Returns
    -------
    tuple
        (X_processed, metadata) where metadata contains processing info
        
    Raises
    ------
    ValueError
        If unknown dimensionality reduction method is specified
    NotImplementedError
        If autoencoder method is requested (not yet implemented)
    """
    metadata = {'original_shape': X.shape}
    
    dim_red_config = config.get('dimensionality_reduction', {})
    if not dim_red_config.get('enabled', False):
        return X, metadata
    
    method = dim_red_config.get('method', 'pca')
    n_components = dim_red_config.get('n_components', 10)
    
    if method == 'pca':
        from sklearn.decomposition import PCA
        actual_components = min(n_components, X.shape[1], X.shape[0])
        pca = PCA(n_components=actual_components)
        X_processed = pca.fit_transform(X)
        metadata['method'] = 'pca'
        metadata['n_components'] = pca.n_components_
        metadata['explained_variance_ratio'] = pca.explained_variance_ratio_.tolist()
        metadata['total_explained_variance'] = float(sum(pca.explained_variance_ratio_))
        return X_processed, metadata
    elif method == 'autoencoder':
        raise NotImplementedError(
            "Autoencoder dimensionality reduction is not yet implemented. "
            "Please use 'pca' method or set 'enabled' to false."
        )
    else:
        raise ValueError(f"Unknown dimensionality reduction method: {method}")


def bootstrap_ate(causal_model, df, graph, n_bootstrap=100, random_state=0):
    """
    Compute ATE with bootstrap confidence intervals.
    
    Parameters
    ----------
    causal_model : gcm.StructuralCausalModel
        Fitted causal model (used as template)
    df : pd.DataFrame
        Original data
    graph : networkx.DiGraph
        Causal graph structure
    n_bootstrap : int
        Number of bootstrap iterations
    random_state : int
        Random seed for reproducibility
        
    Returns
    -------
    dict
        ATE estimates with confidence intervals
    """
    rng = np.random.default_rng(random_state)
    n = len(df)
    ate_samples = []
    
    for i in range(n_bootstrap):
        indices = rng.choice(n, size=n, replace=True)
        df_boot = df.iloc[indices].reset_index(drop=True)
        
        model_boot = gcm.StructuralCausalModel(graph)
        gcm.auto.assign_causal_mechanisms(model_boot, df_boot)
        gcm.fit(model_boot, df_boot)
        
        boot_seed = random_state + i + 1
        
        y1 = _interventional_samples_compat(
            model_boot, {"T": lambda t: 1},
            num_samples_to_draw=n, random_state=boot_seed
        )["Y"]
        if hasattr(y1, 'to_numpy'):
            y1 = y1.to_numpy()
        
        y0 = _interventional_samples_compat(
            model_boot, {"T": lambda t: 0},
            num_samples_to_draw=n, random_state=boot_seed + 1000
        )["Y"]
        if hasattr(y0, 'to_numpy'):
            y0 = y0.to_numpy()
        
        ate_samples.append(float(np.mean(y1 - y0)))
    
    ate_samples = np.array(ate_samples)
    return {
        'mean': float(np.mean(ate_samples)),
        'std': float(np.std(ate_samples)),
        'ci_lower': float(np.percentile(ate_samples, 2.5)),
        'ci_upper': float(np.percentile(ate_samples, 97.5)),
        'samples': ate_samples.tolist()
    }


def summarize_ite_distribution(ite):
    """
    Summarize ITE distribution with quantiles and statistics.
    
    Parameters
    ----------
    ite : numpy.ndarray
        Individual treatment effects
        
    Returns
    -------
    dict
        Summary statistics including mean, std, quantiles, min, max
    """
    ite = np.asarray(ite)
    return {
        'mean': float(np.mean(ite)),
        'std': float(np.std(ite)),
        'quantile_05': float(np.percentile(ite, 5)),
        'quantile_50': float(np.percentile(ite, 50)),
        'quantile_95': float(np.percentile(ite, 95)),
        'min': float(np.min(ite)),
        'max': float(np.max(ite))
    }


def log_run(metadata, output_path='metadata.json'):
    """
    Save run metadata to JSON file.
    
    Parameters
    ----------
    metadata : dict
        Metadata dictionary to save
    output_path : str
        Path to save metadata JSON file
    """
    metadata = metadata.copy()
    metadata['timestamp'] = datetime.now().isoformat()
    
    try:
        import dowhy
        dowhy_version = getattr(dowhy, '__version__', 'unknown')
    except Exception:
        dowhy_version = 'unknown'
    
    try:
        import sklearn
        sklearn_version = sklearn.__version__
    except Exception:
        sklearn_version = 'unknown'
    
    metadata['package_versions'] = {
        'dowhy': dowhy_version,
        'numpy': np.__version__,
        'pandas': pd.__version__,
        'scikit-learn': sklearn_version,
        'networkx': nx.__version__
    }
    
    with open(output_path, 'w') as f:
        json.dump(metadata, f, indent=2, default=str)


def run_dowhy_gcm(X, T, Y, edges_list=None, random_state=0, 
                  config_path=None, dag_path=None, return_metadata=False, **kwargs):
    """
    Run DoWhy GCM causal inference with configurable options.
    
    This function estimates causal effects using DoWhy's Graphical Causal Models.
    It supports flexible DAG specification, high-dimensional data processing,
    and provides uncertainty estimates via bootstrap.
    
    Parameters
    ----------
    X : numpy.ndarray
        Covariate matrix (n, p). Can be high-dimensional; use config for
        dimensionality reduction.
    T : numpy.ndarray
        Treatment vector (n,). Must be binary (0/1).
    Y : numpy.ndarray
        Outcome vector (n,).
    edges_list : list, optional
        List of [parent, child] edges for the DAG. If provided, takes precedence
        over dag_path. This is the primary interface from R.
    random_state : int, default=0
        Random seed for reproducibility.
    config_path : str, optional
        Path to config JSON file. If None, uses default configuration.
    dag_path : str, optional
        Path to DAG JSON file. Only used if edges_list is None.
    return_metadata : bool, default=False
        Whether to include metadata in the return value.
    **kwargs : dict
        Additional arguments:
        - variable_names: list of column names for X
        - bootstrap: bool, whether to run bootstrap ATE estimation
        - n_bootstrap: int, number of bootstrap iterations
        
    Returns
    -------
    dict
        Results dictionary with keys:
        - "ite": Individual treatment effects (numpy array)
        - "ate": Average treatment effect (float)
        - "y0_hat": Counterfactual outcomes under T=0 (numpy array)
        - "y1_hat": Counterfactual outcomes under T=1 (numpy array)
        - "summary": (optional) ATE CI, ITE statistics
        - "samples": (optional) y0, y1, ite arrays as lists
        - "metadata": (optional) input hash, DAG structure, settings
    """
    # Load configuration
    config = load_config(config_path)
    if random_state is not None:
        config['random_seed'] = random_state
    random_state = config.get('random_seed', 0)
    
    # Validate inputs
    X, T, Y = validate_inputs(X, T, Y)
    
    # Process X with optional dimensionality reduction
    X_processed, dim_red_meta = process_X(X, config)
    
    n, p = X_processed.shape
    
    # Handle variable names for X columns
    variable_names = kwargs.get("variable_names", None)
    if variable_names is not None:
        if len(variable_names) != p:
            # If PCA was applied, variable_names won't match
            # Fall back to generic names
            columns = [f"X{i}" for i in range(p)]
        else:
            columns = list(variable_names)
    else:
        columns = [f"X{i}" for i in range(p)]
    
    # Construct DataFrame
    df = pd.DataFrame(X_processed, columns=columns)
    df["T"] = T
    df["Y"] = Y
    
    # Build causal graph
    if edges_list is not None:
        # edges_list from R: list of [parent, child] pairs
        causal_graph = nx.DiGraph()
        causal_graph.add_edges_from(edges_list)
    elif dag_path:
        causal_graph = load_dag(dag_path)
    else:
        # Default DAG structure for backward compatibility
        edges = [("X0", "T"), ("T", "Y")]
        if p > 1:
            edges.append(("X1", "Y"))
        causal_graph = nx.DiGraph(edges)
    
    # Fit structural causal model
    try:
        scm = gcm.StructuralCausalModel(causal_graph)
        gcm.auto.assign_causal_mechanisms(scm, df)
        gcm.fit(scm, df)
    except Exception as e:
        raise RuntimeError(f"Failed to fit causal model: {str(e)}")
    
    # Generate interventional samples (version-compatible)
    n_samples = n
    
    y1_samples = _interventional_samples_compat(
        scm, {"T": lambda t: 1}, num_samples_to_draw=n_samples, random_state=random_state
    )["Y"]
    if hasattr(y1_samples, 'to_numpy'):
        y1 = y1_samples.to_numpy()
    else:
        y1 = np.array(y1_samples)
    
    y0_samples = _interventional_samples_compat(
        scm, {"T": lambda t: 0}, num_samples_to_draw=n_samples, random_state=random_state + 1
    )["Y"]
    if hasattr(y0_samples, 'to_numpy'):
        y0 = y0_samples.to_numpy()
    else:
        y0 = np.array(y0_samples)
    
    # Compute treatment effects
    ite = y1 - y0
    ate = float(np.mean(ite))
    
    # Build result dictionary (backward compatible)
    result = {
        "ate": ate,
        "ite": ite,
        "y1_hat": y1,
        "y0_hat": y0,
    }
    
    # Optional bootstrap ATE estimation
    do_bootstrap = kwargs.get('bootstrap', config.get('bootstrap', {}).get('enabled', False))
    n_bootstrap = kwargs.get('n_bootstrap', config.get('bootstrap', {}).get('n_iterations', 100))
    
    if do_bootstrap:
        ate_bootstrap = bootstrap_ate(scm, df, causal_graph, n_bootstrap, random_state)
        ite_summary = summarize_ite_distribution(ite)
        
        result["summary"] = {
            "ate": ate,
            "ate_ci_lower": ate_bootstrap['ci_lower'],
            "ate_ci_upper": ate_bootstrap['ci_upper'],
            "ate_std": ate_bootstrap['std'],
            "ite_mean": ite_summary['mean'],
            "ite_std": ite_summary['std'],
            "ite_quantiles": {
                "q05": ite_summary['quantile_05'],
                "q50": ite_summary['quantile_50'],
                "q95": ite_summary['quantile_95']
            }
        }
        result["samples"] = {
            "y0": y0.tolist(),
            "y1": y1.tolist(),
            "ite": ite.tolist(),
            "ate_bootstrap": ate_bootstrap['samples']
        }
    else:
        # Basic summary without bootstrap
        ite_summary = summarize_ite_distribution(ite)
        result["summary"] = {
            "ate": ate,
            "ite_mean": ite_summary['mean'],
            "ite_std": ite_summary['std'],
            "ite_quantiles": {
                "q05": ite_summary['quantile_05'],
                "q50": ite_summary['quantile_50'],
                "q95": ite_summary['quantile_95']
            }
        }
        result["samples"] = {
            "y0": y0.tolist(),
            "y1": y1.tolist(),
            "ite": ite.tolist()
        }
    
    # Optional metadata
    if return_metadata or config.get('output_metadata', False):
        input_data = np.concatenate([X.flatten(), T.flatten(), Y.flatten()])
        input_hash = hashlib.sha256(input_data.tobytes()).hexdigest()
        
        result["metadata"] = {
            "input_hash": input_hash,
            "dag_structure": {
                "nodes": list(causal_graph.nodes()),
                "edges": list(causal_graph.edges())
            },
            "intervention_conditions": ["do(T=0)", "do(T=1)"],
            "random_seed": random_state,
            "n_samples": n,
            "n_features_original": X.shape[1],
            "n_features_processed": p,
            "dimensionality_reduction": dim_red_meta,
            "bootstrap_enabled": do_bootstrap
        }
        
        if do_bootstrap:
            result["metadata"]["bootstrap_iterations"] = n_bootstrap
        
        if config.get('output_metadata', False):
            output_path = config.get('metadata_output_path', 'metadata.json')
            log_run(result["metadata"], output_path)
    
    return result


if __name__ == "__main__":
    pass

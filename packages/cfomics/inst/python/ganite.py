# inst/python/ganite.py
# GANITE implementation for cfomics package
#
# Reference: Jinsung Yoon, James Jordon, Mihaela van der Schaar,
# "GANITE: Estimation of Individualized Treatment Effects using Generative Adversarial Nets",
# International Conference on Learning Representations (ICLR), 2018.
# Paper link: https://openreview.net/forum?id=ByKWUeWA-

import numpy as np
import pandas as pd


def _batch_generator(x, t, y, batch_size, rng):
    """Generate mini-batch with x, t, and y.
    
    Args:
        x: features (n, p)
        t: treatments (n,)
        y: observed labels (n,)
        batch_size: mini batch size
        rng: numpy random generator
        
    Returns:
        X_mb, T_mb, Y_mb: mini-batch data
    """
    n = x.shape[0]
    batch_idx = rng.choice(n, size=min(batch_size, n), replace=False)
    
    X_mb = x[batch_idx, :]
    T_mb = t[batch_idx].reshape(-1, 1)
    Y_mb = y[batch_idx].reshape(-1, 1)
    return X_mb, T_mb, Y_mb


class GANITEModel:
    """GANITE model using TensorFlow 2.x/Keras.
    
    This implements the GANITE algorithm for Individual Treatment Effect estimation
    using Generative Adversarial Networks.
    """
    
    def __init__(self, input_dim, h_dim=100, alpha=1.0, random_state=0):
        """Initialize GANITE model.
        
        Args:
            input_dim: dimension of input features
            h_dim: hidden layer dimension
            alpha: hyperparameter to balance GAN loss and factual loss
            random_state: random seed for reproducibility
        """
        self.input_dim = input_dim
        self.h_dim = h_dim
        self.alpha = alpha
        self.random_state = random_state
        
        try:
            import tensorflow as tf
            tf.random.set_seed(random_state)
            self.tf = tf
        except ImportError:
            raise ImportError(
                "TensorFlow is required for GANITE. "
                "Install it with: pip install tensorflow"
            )
        
        self._build_networks()
    
    def _build_networks(self):
        """Build Generator, Discriminator, and Inference networks."""
        tf = self.tf
        from tensorflow import keras
        from tensorflow.keras import layers
        
        # Generator: takes X, T, Y and outputs potential outcomes (Y0, Y1)
        self.generator = keras.Sequential([
            layers.Dense(self.h_dim, activation='relu', 
                        input_shape=(self.input_dim + 2,),
                        kernel_initializer='glorot_uniform'),
            layers.Dense(self.h_dim, activation='relu',
                        kernel_initializer='glorot_uniform'),
        ], name='generator_base')
        
        # Generator output heads
        self.gen_head_y0 = keras.Sequential([
            layers.Dense(self.h_dim, activation='relu',
                        kernel_initializer='glorot_uniform'),
            layers.Dense(1, kernel_initializer='glorot_uniform')
        ], name='gen_head_y0')
        
        self.gen_head_y1 = keras.Sequential([
            layers.Dense(self.h_dim, activation='relu',
                        kernel_initializer='glorot_uniform'),
            layers.Dense(1, kernel_initializer='glorot_uniform')
        ], name='gen_head_y1')
        
        # Discriminator: takes X, Y0, Y1 and predicts treatment
        self.discriminator = keras.Sequential([
            layers.Dense(self.h_dim, activation='relu',
                        input_shape=(self.input_dim + 2,),
                        kernel_initializer='glorot_uniform'),
            layers.Dense(self.h_dim, activation='relu',
                        kernel_initializer='glorot_uniform'),
            layers.Dense(1, kernel_initializer='glorot_uniform')
        ], name='discriminator')
        
        # Inference network: takes X and outputs potential outcomes (Y0, Y1)
        self.inference = keras.Sequential([
            layers.Dense(self.h_dim, activation='relu',
                        input_shape=(self.input_dim,),
                        kernel_initializer='glorot_uniform'),
            layers.Dense(self.h_dim, activation='relu',
                        kernel_initializer='glorot_uniform'),
        ], name='inference_base')
        
        # Inference output heads
        self.inf_head_y0 = keras.Sequential([
            layers.Dense(self.h_dim, activation='relu',
                        kernel_initializer='glorot_uniform'),
            layers.Dense(1, kernel_initializer='glorot_uniform')
        ], name='inf_head_y0')
        
        self.inf_head_y1 = keras.Sequential([
            layers.Dense(self.h_dim, activation='relu',
                        kernel_initializer='glorot_uniform'),
            layers.Dense(1, kernel_initializer='glorot_uniform')
        ], name='inf_head_y1')
        
        # Optimizers
        self.g_optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)
        self.d_optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)
        self.i_optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)
    
    def _generator_forward(self, x, t, y):
        """Forward pass through generator."""
        tf = self.tf
        inputs = tf.concat([x, t, y], axis=1)
        h = self.generator(inputs)
        y0_logit = self.gen_head_y0(h)
        y1_logit = self.gen_head_y1(h)
        return y0_logit, y1_logit
    
    def _inference_forward(self, x):
        """Forward pass through inference network."""
        h = self.inference(x)
        y0_logit = self.inf_head_y0(h)
        y1_logit = self.inf_head_y1(h)
        return y0_logit, y1_logit
    
    def _discriminator_forward(self, x, t, y, y_tilde):
        """Forward pass through discriminator."""
        tf = self.tf
        # Combine factual and counterfactual outcomes
        y0_tilde, y1_tilde = y_tilde[:, 0:1], y_tilde[:, 1:2]
        input0 = (1.0 - t) * y + t * y0_tilde  # Y0: factual if T=0, counterfactual if T=1
        input1 = t * y + (1.0 - t) * y1_tilde  # Y1: factual if T=1, counterfactual if T=0
        inputs = tf.concat([x, input0, input1], axis=1)
        return self.discriminator(inputs)
    
    @property
    def generator_variables(self):
        return (self.generator.trainable_variables + 
                self.gen_head_y0.trainable_variables + 
                self.gen_head_y1.trainable_variables)
    
    @property
    def inference_variables(self):
        return (self.inference.trainable_variables + 
                self.inf_head_y0.trainable_variables + 
                self.inf_head_y1.trainable_variables)
    
    def train_step_gd(self, x, t, y):
        """Single training step for Generator and Discriminator."""
        tf = self.tf
        
        # Train Discriminator
        with tf.GradientTape() as tape:
            y0_logit, y1_logit = self._generator_forward(x, t, y)
            y_tilde = tf.concat([tf.sigmoid(y0_logit), tf.sigmoid(y1_logit)], axis=1)
            d_logit = self._discriminator_forward(x, t, y, y_tilde)
            d_loss = tf.reduce_mean(
                tf.nn.sigmoid_cross_entropy_with_logits(labels=t, logits=d_logit)
            )
        d_grads = tape.gradient(d_loss, self.discriminator.trainable_variables)
        self.d_optimizer.apply_gradients(
            zip(d_grads, self.discriminator.trainable_variables)
        )
        
        # Train Generator
        with tf.GradientTape() as tape:
            y0_logit, y1_logit = self._generator_forward(x, t, y)
            y_tilde = tf.concat([tf.sigmoid(y0_logit), tf.sigmoid(y1_logit)], axis=1)
            d_logit = self._discriminator_forward(x, t, y, y_tilde)
            
            # GAN loss (fool discriminator)
            g_loss_gan = -tf.reduce_mean(
                tf.nn.sigmoid_cross_entropy_with_logits(labels=t, logits=d_logit)
            )
            
            # Factual loss
            factual_logit = t * y1_logit + (1.0 - t) * y0_logit
            g_loss_factual = tf.reduce_mean(
                tf.nn.sigmoid_cross_entropy_with_logits(labels=y, logits=factual_logit)
            )
            
            g_loss = g_loss_factual + self.alpha * g_loss_gan
        
        g_grads = tape.gradient(g_loss, self.generator_variables)
        self.g_optimizer.apply_gradients(zip(g_grads, self.generator_variables))
        
        return float(d_loss.numpy()), float(g_loss.numpy())
    
    def train_step_inference(self, x, t, y):
        """Single training step for Inference network."""
        tf = self.tf
        
        with tf.GradientTape() as tape:
            # Get generator outputs (for counterfactual labels)
            y0_logit_g, y1_logit_g = self._generator_forward(x, t, y)
            y0_tilde = tf.sigmoid(y0_logit_g)
            y1_tilde = tf.sigmoid(y1_logit_g)
            
            # Get inference outputs
            y0_logit_i, y1_logit_i = self._inference_forward(x)
            
            # Labels for inference network
            # For Y1: use factual Y if T=1, use generator's Y1 if T=0
            label_y1 = t * y + (1.0 - t) * y1_tilde
            # For Y0: use factual Y if T=0, use generator's Y0 if T=1
            label_y0 = (1.0 - t) * y + t * y0_tilde
            
            i_loss_y1 = tf.reduce_mean(
                tf.nn.sigmoid_cross_entropy_with_logits(labels=label_y1, logits=y1_logit_i)
            )
            i_loss_y0 = tf.reduce_mean(
                tf.nn.sigmoid_cross_entropy_with_logits(labels=label_y0, logits=y0_logit_i)
            )
            i_loss = i_loss_y0 + i_loss_y1
        
        i_grads = tape.gradient(i_loss, self.inference_variables)
        self.i_optimizer.apply_gradients(zip(i_grads, self.inference_variables))
        
        return float(i_loss.numpy())
    
    def fit(self, X, T, Y, iterations=10000, batch_size=256, verbose=True):
        """Train the GANITE model.
        
        Args:
            X: features (n, p)
            T: treatments (n,)
            Y: outcomes (n,)
            iterations: number of training iterations
            batch_size: mini-batch size
            verbose: whether to print progress
        """
        tf = self.tf
        rng = np.random.default_rng(self.random_state)
        
        # Normalize Y to [0, 1] for sigmoid outputs
        self.y_min = Y.min()
        self.y_max = Y.max()
        self.y_range = self.y_max - self.y_min
        if self.y_range < 1e-8:
            self.y_range = 1.0
        Y_norm = (Y - self.y_min) / self.y_range
        
        # Phase 1: Train Generator and Discriminator
        if verbose:
            print('Training Generator and Discriminator...')
        
        for it in range(iterations):
            # Train discriminator twice per generator update
            for _ in range(2):
                X_mb, T_mb, Y_mb = _batch_generator(X, T, Y_norm, batch_size, rng)
                X_mb = tf.cast(X_mb, tf.float32)
                T_mb = tf.cast(T_mb, tf.float32)
                Y_mb = tf.cast(Y_mb, tf.float32)
                d_loss, _ = self.train_step_gd(X_mb, T_mb, Y_mb)
            
            X_mb, T_mb, Y_mb = _batch_generator(X, T, Y_norm, batch_size, rng)
            X_mb = tf.cast(X_mb, tf.float32)
            T_mb = tf.cast(T_mb, tf.float32)
            Y_mb = tf.cast(Y_mb, tf.float32)
            d_loss, g_loss = self.train_step_gd(X_mb, T_mb, Y_mb)
            
            if verbose and it % 1000 == 0:
                print(f'Iteration {it}/{iterations}, D loss: {d_loss:.4f}, G loss: {g_loss:.4f}')
        
        # Phase 2: Train Inference network
        if verbose:
            print('Training Inference network...')
        
        for it in range(iterations):
            X_mb, T_mb, Y_mb = _batch_generator(X, T, Y_norm, batch_size, rng)
            X_mb = tf.cast(X_mb, tf.float32)
            T_mb = tf.cast(T_mb, tf.float32)
            Y_mb = tf.cast(Y_mb, tf.float32)
            i_loss = self.train_step_inference(X_mb, T_mb, Y_mb)
            
            if verbose and it % 1000 == 0:
                print(f'Iteration {it}/{iterations}, I loss: {i_loss:.4f}')
    
    def predict(self, X):
        """Predict potential outcomes for given features.
        
        Args:
            X: features (n, p)
            
        Returns:
            y_hat: predicted potential outcomes (n, 2) where [:, 0] is Y0 and [:, 1] is Y1
        """
        tf = self.tf
        X_tf = tf.cast(X, tf.float32)
        y0_logit, y1_logit = self._inference_forward(X_tf)
        y0_hat = tf.sigmoid(y0_logit).numpy()
        y1_hat = tf.sigmoid(y1_logit).numpy()
        
        # Denormalize
        y0_hat = y0_hat * self.y_range + self.y_min
        y1_hat = y1_hat * self.y_range + self.y_min
        
        return np.concatenate([y0_hat, y1_hat], axis=1)


def run_ganite(X, T, Y, edges_list=None, variable_names=None, random_state=0, **kwargs):
    """
    GANITE implementation for individual treatment effect estimation.
    
    Parameters
    ----------
    X : numpy array, shape (n, p)
        Covariate matrix
    T : numpy array, shape (n,)
        Treatment vector (binary: 0 or 1)
    Y : numpy array, shape (n,)
        Outcome vector
    edges_list : list of [parent, child], optional
        DAG edges (not directly used by GANITE but kept for API consistency)
    variable_names : list of str, optional
        Names of covariates in X
    random_state : int
        Random seed for reproducibility
    **kwargs : dict
        Additional parameters for GANITE:
        - h_dim: hidden layer dimension (default: 100)
        - alpha: GAN loss weight (default: 1.0)
        - iterations: training iterations (default: 10000)
        - batch_size: mini-batch size (default: 256)
        - verbose: print progress (default: False)
    
    Returns
    -------
    dict with keys:
        - "ate": float, average treatment effect
        - "ite": numpy array, individual treatment effects
        - "y1_hat": numpy array, predicted outcomes under treatment
        - "y0_hat": numpy array, predicted outcomes under control
        - "summary": dict with summary statistics
        - "samples": dict with sample arrays
    """
    # Convert inputs to numpy arrays
    if not isinstance(X, np.ndarray):
        X = np.asarray(X)
    if not isinstance(T, np.ndarray):
        T = np.asarray(T)
    if not isinstance(Y, np.ndarray):
        Y = np.asarray(Y)
    
    # Ensure correct shapes
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    T = T.flatten()
    Y = Y.flatten()
    
    # Validate inputs
    n, p = X.shape
    if T.shape[0] != n:
        raise ValueError(f"Shape mismatch: X has {n} samples but T has {T.shape[0]}")
    if Y.shape[0] != n:
        raise ValueError(f"Shape mismatch: X has {n} samples but Y has {Y.shape[0]}")
    
    unique_T = np.unique(T)
    if not np.all(np.isin(unique_T, [0, 1])):
        raise ValueError(f"T must be binary (0/1), got unique values: {unique_T}")
    
    # Extract GANITE parameters
    h_dim = kwargs.get('h_dim', 100)
    alpha = kwargs.get('alpha', 1.0)
    iterations = kwargs.get('iterations', 10000)
    batch_size = kwargs.get('batch_size', 256)
    verbose = kwargs.get('verbose', False)
    
    # Build and train model
    model = GANITEModel(
        input_dim=p,
        h_dim=h_dim,
        alpha=alpha,
        random_state=random_state
    )
    
    model.fit(X, T, Y, iterations=iterations, batch_size=batch_size, verbose=verbose)
    
    # Predict potential outcomes
    y_hat = model.predict(X)
    y0_hat = y_hat[:, 0]
    y1_hat = y_hat[:, 1]
    
    # Compute treatment effects
    ite = y1_hat - y0_hat
    ate = float(np.mean(ite))
    
    # Summary statistics
    summary = {
        'ate': ate,
        'ite_mean': float(np.mean(ite)),
        'ite_std': float(np.std(ite)),
        'ite_quantiles': {
            'q05': float(np.percentile(ite, 5)),
            'q50': float(np.percentile(ite, 50)),
            'q95': float(np.percentile(ite, 95))
        }
    }
    
    # Samples
    samples = {
        'y0': y0_hat.tolist(),
        'y1': y1_hat.tolist(),
        'ite': ite.tolist()
    }
    
    return {
        'ate': ate,
        'ite': ite,
        'y1_hat': y1_hat,
        'y0_hat': y0_hat,
        'summary': summary,
        'samples': samples
    }

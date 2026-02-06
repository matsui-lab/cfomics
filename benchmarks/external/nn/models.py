# benchmarks/external/nn/models.py
"""Neural network models for causal inference: TARNet, CFRNet, DragonNet

These models are for PAPER COMPARISON ONLY (not part of cfomics package).
They implement standard causal inference NN architectures from the literature:

- TARNet: Two-head architecture for Y(0) and Y(1)
  Reference: Shalit et al. (2017) "Estimating individual treatment effect"

- CFRNet: TARNet + IPM regularization for balanced representations
  Reference: Shalit et al. (2017) "Estimating individual treatment effect"

- DragonNet: TARNet + propensity score co-learning
  Reference: Shi et al. (2019) "Adapting Neural Networks for the Estimation
  of Treatment Effects"
"""

import torch
import torch.nn as nn
import numpy as np
from typing import Tuple, Optional, Dict, Any


class RepresentationNet(nn.Module):
    """Shared representation network.

    Maps input covariates to a learned representation space.
    Uses ELU activation for smooth gradients.

    Args:
        input_dim: Dimension of input covariates
        hidden_dim: Dimension of hidden layers (default: 200)
        n_layers: Number of hidden layers (default: 3)
    """

    def __init__(self, input_dim: int, hidden_dim: int = 200, n_layers: int = 3):
        super().__init__()
        layers = [nn.Linear(input_dim, hidden_dim), nn.ELU()]
        for _ in range(n_layers - 1):
            layers.extend([nn.Linear(hidden_dim, hidden_dim), nn.ELU()])
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class TARNet(nn.Module):
    """Treatment-Agnostic Representation Network.

    Two-head architecture that learns a shared representation and
    separate outcome heads for treated and control groups.

    Reference:
        Shalit, Johansson, Sontag (2017). "Estimating individual treatment
        effect: generalization bounds and algorithms." ICML.

    Args:
        input_dim: Dimension of input covariates
        hidden_dim: Dimension of hidden layers (default: 200)
        n_layers: Number of representation layers (default: 3)
    """

    def __init__(self, input_dim: int, hidden_dim: int = 200, n_layers: int = 3):
        super().__init__()
        self.repr = RepresentationNet(input_dim, hidden_dim, n_layers)
        self.head0 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim), nn.ELU(),
            nn.Linear(hidden_dim, 1)
        )
        self.head1 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim), nn.ELU(),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, x: torch.Tensor, t: torch.Tensor) -> Tuple[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """Forward pass.

        Args:
            x: Covariate tensor of shape (batch_size, input_dim)
            t: Treatment tensor of shape (batch_size,) with values 0 or 1

        Returns:
            Tuple of (y_pred, (y0, y1)) where:
                - y_pred: Factual prediction based on actual treatment
                - y0: Counterfactual prediction under control
                - y1: Counterfactual prediction under treatment
        """
        phi = self.repr(x)
        y0 = self.head0(phi).squeeze(-1)
        y1 = self.head1(phi).squeeze(-1)
        y_pred = t * y1 + (1 - t) * y0
        return y_pred, (y0, y1)


class CFRNet(TARNet):
    """Counterfactual Regression Network with IPM regularization.

    Extends TARNet with Integral Probability Metric (IPM) regularization
    to encourage balanced representations between treatment groups.
    Uses Wasserstein-1 approximation via mean matching.

    Reference:
        Shalit, Johansson, Sontag (2017). "Estimating individual treatment
        effect: generalization bounds and algorithms." ICML.

    Args:
        input_dim: Dimension of input covariates
        hidden_dim: Dimension of hidden layers (default: 200)
        n_layers: Number of representation layers (default: 3)
        alpha: Weight for IPM regularization term (default: 1.0)
    """

    def __init__(self, input_dim: int, hidden_dim: int = 200, n_layers: int = 3, alpha: float = 1.0):
        super().__init__(input_dim, hidden_dim, n_layers)
        self.alpha = alpha

    def ipm_loss(self, phi: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        """Wasserstein-1 (MMD) approximation for IPM.

        Computes the L2 distance between mean representations of
        treated and control groups as a proxy for the Wasserstein distance.

        Args:
            phi: Representation tensor of shape (batch_size, hidden_dim)
            t: Treatment tensor of shape (batch_size,)

        Returns:
            Scalar IPM loss value
        """
        t = t.bool()
        phi0 = phi[~t]
        phi1 = phi[t]
        if len(phi0) == 0 or len(phi1) == 0:
            return torch.tensor(0.0, device=phi.device)
        mean0 = phi0.mean(dim=0)
        mean1 = phi1.mean(dim=0)
        return torch.norm(mean0 - mean1, p=2)


class DragonNet(nn.Module):
    """DragonNet: propensity score co-learning (Simplified Version).

    Based on Shi et al. (2019) "Adapting Neural Networks for Treatment Effects".

    NOTE: This is a simplified implementation that includes joint propensity
    learning but omits the targeted regularization loss from the original paper.
    For full DragonNet, see the official implementation at:
    https://github.com/claudiashi57/dragonnet

    The simplified version still provides benefits from representation sharing
    between outcome and propensity models, but may have slightly higher bias
    than the full targeted regularization approach.

    Args:
        input_dim: Number of input features
        hidden_dim: Hidden layer dimension (default: 200)
        n_layers: Number of representation layers (default: 3)
    """

    def __init__(self, input_dim: int, hidden_dim: int = 200, n_layers: int = 3):
        super().__init__()
        self.repr = RepresentationNet(input_dim, hidden_dim, n_layers)
        self.head0 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim), nn.ELU(),
            nn.Linear(hidden_dim, 1)
        )
        self.head1 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim), nn.ELU(),
            nn.Linear(hidden_dim, 1)
        )
        self.propensity_head = nn.Sequential(
            nn.Linear(hidden_dim, 1),
            nn.Sigmoid()
        )

    def forward(self, x: torch.Tensor, t: torch.Tensor) -> Tuple[torch.Tensor, Tuple[torch.Tensor, torch.Tensor], torch.Tensor]:
        """Forward pass.

        Args:
            x: Covariate tensor of shape (batch_size, input_dim)
            t: Treatment tensor of shape (batch_size,) with values 0 or 1

        Returns:
            Tuple of (y_pred, (y0, y1), ps) where:
                - y_pred: Factual prediction based on actual treatment
                - y0: Counterfactual prediction under control
                - y1: Counterfactual prediction under treatment
                - ps: Estimated propensity score P(T=1|X)
        """
        phi = self.repr(x)
        y0 = self.head0(phi).squeeze(-1)
        y1 = self.head1(phi).squeeze(-1)
        ps = self.propensity_head(phi).squeeze(-1)
        y_pred = t * y1 + (1 - t) * y0
        return y_pred, (y0, y1), ps


def train_model(
    model: nn.Module,
    X: np.ndarray,
    T: np.ndarray,
    Y: np.ndarray,
    n_epochs: int = 300,
    lr: float = 1e-3,
    batch_size: int = 64,
    alpha: float = 1.0,
    verbose: bool = False
) -> nn.Module:
    """Train a causal NN model.

    Args:
        model: TARNet, CFRNet, or DragonNet instance
        X: Covariate array of shape (n, p)
        T: Treatment array of shape (n,) with values 0 or 1
        Y: Outcome array of shape (n,)
        n_epochs: Number of training epochs (default: 300)
        lr: Learning rate (default: 1e-3)
        batch_size: Mini-batch size (default: 64)
        alpha: IPM regularization weight for CFRNet (default: 1.0)
        verbose: Print training progress (default: False)

    Returns:
        Trained model
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)

    X_t = torch.tensor(X, dtype=torch.float32, device=device)
    T_t = torch.tensor(T, dtype=torch.float32, device=device)
    Y_t = torch.tensor(Y, dtype=torch.float32, device=device)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    n = len(Y)

    for epoch in range(n_epochs):
        model.train()
        perm = torch.randperm(n)
        epoch_loss = 0.0
        n_batches = 0

        for i in range(0, n, batch_size):
            idx = perm[i:i+batch_size]
            x_b, t_b, y_b = X_t[idx], T_t[idx], Y_t[idx]

            optimizer.zero_grad()

            if isinstance(model, DragonNet):
                # Simplified DragonNet loss: outcome + propensity (no targeted regularization)
                y_pred, (y0, y1), ps = model(x_b, t_b)
                loss_y = nn.MSELoss()(y_pred, y_b)
                loss_ps = nn.BCELoss()(ps, t_b)
                loss = loss_y + loss_ps
            elif isinstance(model, CFRNet):
                y_pred, (y0, y1) = model(x_b, t_b)
                loss_y = nn.MSELoss()(y_pred, y_b)
                phi = model.repr(x_b)
                loss_ipm = model.ipm_loss(phi, t_b)
                loss = loss_y + alpha * loss_ipm
            else:  # TARNet
                y_pred, (y0, y1) = model(x_b, t_b)
                loss = nn.MSELoss()(y_pred, y_b)

            loss.backward()
            optimizer.step()

            epoch_loss += loss.item()
            n_batches += 1

        if verbose and (epoch + 1) % 50 == 0:
            avg_loss = epoch_loss / n_batches
            print(f"Epoch {epoch + 1}/{n_epochs}, Loss: {avg_loss:.4f}")

    return model


def predict_ite(model: nn.Module, X: np.ndarray) -> Dict[str, Any]:
    """Predict ITE using trained model.

    Args:
        model: Trained TARNet, CFRNet, or DragonNet
        X: Covariate array of shape (n, p)

    Returns:
        Dictionary with:
            - ite: Individual treatment effects array
            - ate: Average treatment effect scalar
            - y0_hat: Predicted outcomes under control
            - y1_hat: Predicted outcomes under treatment
    """
    device = next(model.parameters()).device
    model.eval()

    with torch.no_grad():
        X_t = torch.tensor(X, dtype=torch.float32, device=device)
        T0 = torch.zeros(len(X), device=device)
        T1 = torch.ones(len(X), device=device)

        if isinstance(model, DragonNet):
            _, (y0, _), _ = model(X_t, T0)
            _, (_, y1), _ = model(X_t, T1)
        else:
            _, (y0, y1) = model(X_t, T0)

        ite = (y1 - y0).cpu().numpy()
        y0_hat = y0.cpu().numpy()
        y1_hat = y1.cpu().numpy()

    return {
        "ite": ite,
        "ate": float(np.mean(ite)),
        "y0_hat": y0_hat,
        "y1_hat": y1_hat
    }


if __name__ == "__main__":
    # Quick test with synthetic data
    np.random.seed(42)
    n, p = 500, 10
    X = np.random.randn(n, p)
    T = (np.random.rand(n) > 0.5).astype(float)
    Y = X[:, 0] + T * (0.5 + X[:, 1]) + np.random.randn(n) * 0.1

    print("Testing TARNet...")
    model = TARNet(p, hidden_dim=100, n_layers=2)
    model = train_model(model, X, T, Y, n_epochs=100, verbose=True)
    result = predict_ite(model, X)
    print(f"TARNet ATE: {result['ate']:.3f}")

    print("\nTesting CFRNet...")
    model = CFRNet(p, hidden_dim=100, n_layers=2, alpha=1.0)
    model = train_model(model, X, T, Y, n_epochs=100, verbose=True)
    result = predict_ite(model, X)
    print(f"CFRNet ATE: {result['ate']:.3f}")

    print("\nTesting DragonNet...")
    model = DragonNet(p, hidden_dim=100, n_layers=2)
    model = train_model(model, X, T, Y, n_epochs=100, verbose=True)
    result = predict_ite(model, X)
    print(f"DragonNet ATE: {result['ate']:.3f}")

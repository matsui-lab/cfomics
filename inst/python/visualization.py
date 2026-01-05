import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_ite_distribution(ite, save_path=None, figsize=(10, 6)):
    """
    Plot ITE density distribution.
    
    Parameters
    ----------
    ite : array-like
        Individual treatment effects
    save_path : str, optional
        Path to save the plot. If None, plot is not saved.
    figsize : tuple, default=(10, 6)
        Figure size (width, height)
        
    Returns
    -------
    tuple
        (fig, ax) matplotlib figure and axes objects
    """
    ite = np.asarray(ite)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    ax.hist(ite, bins=50, density=True, alpha=0.7, edgecolor='black', color='steelblue')
    
    mean_ite = np.mean(ite)
    median_ite = np.median(ite)
    
    ax.axvline(mean_ite, color='red', linestyle='--', linewidth=2,
               label=f'Mean ITE: {mean_ite:.3f}')
    ax.axvline(median_ite, color='orange', linestyle=':', linewidth=2,
               label=f'Median ITE: {median_ite:.3f}')
    
    q05 = np.percentile(ite, 5)
    q95 = np.percentile(ite, 95)
    ax.axvspan(q05, q95, alpha=0.2, color='green', label=f'90% CI: [{q05:.3f}, {q95:.3f}]')
    
    ax.set_xlabel('Individual Treatment Effect (ITE)', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Distribution of Individual Treatment Effects', fontsize=14)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig, ax


def plot_outcome_shift(y0, y1, save_path=None, figsize=(14, 6)):
    """
    Plot comparison of Y under T=0 vs T=1.
    
    Parameters
    ----------
    y0 : array-like
        Outcomes under control (T=0)
    y1 : array-like
        Outcomes under treatment (T=1)
    save_path : str, optional
        Path to save the plot. If None, plot is not saved.
    figsize : tuple, default=(14, 6)
        Figure size (width, height)
        
    Returns
    -------
    tuple
        (fig, axes) matplotlib figure and axes objects
    """
    y0 = np.asarray(y0)
    y1 = np.asarray(y1)
    
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    axes[0].hist(y0, bins=50, alpha=0.6, label=f'T=0 (mean={np.mean(y0):.3f})',
                 edgecolor='black', color='blue')
    axes[0].hist(y1, bins=50, alpha=0.6, label=f'T=1 (mean={np.mean(y1):.3f})',
                 edgecolor='black', color='red')
    axes[0].set_xlabel('Outcome Y', fontsize=12)
    axes[0].set_ylabel('Frequency', fontsize=12)
    axes[0].set_title('Outcome Distribution: T=0 vs T=1', fontsize=14)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].scatter(y0, y1, alpha=0.5, s=10, c='steelblue')
    min_val = min(y0.min(), y1.min())
    max_val = max(y0.max(), y1.max())
    axes[1].plot([min_val, max_val], [min_val, max_val],
                 'r--', linewidth=2, label='y1=y0 (no effect)')
    axes[1].set_xlabel('Y under T=0', fontsize=12)
    axes[1].set_ylabel('Y under T=1', fontsize=12)
    axes[1].set_title('Individual Outcome Shift', fontsize=14)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig, axes


def plot_uplift_curve(ite, save_path=None, figsize=(10, 6)):
    """
    Plot uplift curve (cumulative ITE sorted by effect size).
    
    The uplift curve shows the cumulative treatment effect when targeting
    individuals in order of their predicted treatment effect (highest first).
    
    Parameters
    ----------
    ite : array-like
        Individual treatment effects
    save_path : str, optional
        Path to save the plot. If None, plot is not saved.
    figsize : tuple, default=(10, 6)
        Figure size (width, height)
        
    Returns
    -------
    tuple
        (fig, ax) matplotlib figure and axes objects
    """
    ite = np.asarray(ite)
    
    ite_sorted = np.sort(ite)[::-1]
    cumulative_ite = np.cumsum(ite_sorted)
    percentiles = np.arange(1, len(ite) + 1) / len(ite) * 100
    
    fig, ax = plt.subplots(figsize=figsize)
    
    ax.plot(percentiles, cumulative_ite, linewidth=2, color='steelblue',
            label='Cumulative Treatment Effect')
    
    ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    
    random_cumulative = np.cumsum(np.full(len(ite), np.mean(ite)))
    ax.plot(percentiles, random_cumulative, 'r--', linewidth=1.5,
            label='Random Targeting', alpha=0.7)
    
    ax.set_xlabel('Population Percentile (%)', fontsize=12)
    ax.set_ylabel('Cumulative Treatment Effect', fontsize=12)
    ax.set_title('Uplift Curve', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig, ax


def plot_ate_bootstrap_distribution(bootstrap_samples, ate_point=None, save_path=None, figsize=(10, 6)):
    """
    Plot the bootstrap distribution of ATE estimates.
    
    Parameters
    ----------
    bootstrap_samples : array-like
        Bootstrap ATE samples
    ate_point : float, optional
        Point estimate of ATE to highlight
    save_path : str, optional
        Path to save the plot. If None, plot is not saved.
    figsize : tuple, default=(10, 6)
        Figure size (width, height)
        
    Returns
    -------
    tuple
        (fig, ax) matplotlib figure and axes objects
    """
    bootstrap_samples = np.asarray(bootstrap_samples)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    ax.hist(bootstrap_samples, bins=30, density=True, alpha=0.7,
            edgecolor='black', color='steelblue')
    
    ci_lower = np.percentile(bootstrap_samples, 2.5)
    ci_upper = np.percentile(bootstrap_samples, 97.5)
    ax.axvspan(ci_lower, ci_upper, alpha=0.2, color='green',
               label=f'95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]')
    
    if ate_point is not None:
        ax.axvline(ate_point, color='red', linestyle='--', linewidth=2,
                   label=f'Point Estimate: {ate_point:.3f}')
    
    ax.axvline(np.mean(bootstrap_samples), color='orange', linestyle=':',
               linewidth=2, label=f'Bootstrap Mean: {np.mean(bootstrap_samples):.3f}')
    
    ax.set_xlabel('Average Treatment Effect (ATE)', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Bootstrap Distribution of ATE', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig, ax

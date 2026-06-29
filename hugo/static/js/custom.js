// Custom JavaScript for Mini Omics Data Portal

document.addEventListener('DOMContentLoaded', function() {
    // Initialize tooltips
    var tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="tooltip"]'));
    var tooltipList = tooltipTriggerList.map(function (tooltipTriggerEl) {
        return new bootstrap.Tooltip(tooltipTriggerEl);
    });

    // Initialize popovers
    var popoverTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="popover"]'));
    var popoverList = popoverTriggerList.map(function (popoverTriggerEl) {
        return new bootstrap.Popover(popoverTriggerEl);
    });

    // Enhanced iframe loading management
    initializeIframeManagement();
    
    // Initialize search functionality
    initializeSearch();

    // Smooth scrolling for anchor links
    document.querySelectorAll('a[href^="#"]').forEach(anchor => {
        anchor.addEventListener('click', function (e) {
            e.preventDefault();
            const target = document.querySelector(this.getAttribute('href'));
            if (target) {
                target.scrollIntoView({
                    behavior: 'smooth',
                    block: 'start'
                });
            }
        });
    });

    // Add loading state to buttons
    document.querySelectorAll('.btn').forEach(button => {
        button.addEventListener('click', function() {
            if (this.classList.contains('btn-primary')) {
                this.innerHTML = '<i class="bi bi-arrow-clockwise spin me-2"></i>Loading...';
                setTimeout(() => {
                    this.innerHTML = '<i class="bi bi-arrow-right-circle me-2"></i>Explore Data';
                }, 1000);
            }
        });
    });

    // Format large numbers for display
    function formatNumber(num) {
        const value = Number(num);

        if (Number.isNaN(value)) {
            return '0';
        }

        if (value < 10000) {
            return Math.round(value).toString();
        }

        let scaled;
        let suffix;

        if (value < 1000000) {
            scaled = value / 1000;
            suffix = 'K';
        } else if (value < 1000000000) {
            scaled = value / 1000000;
            suffix = 'M';
        } else if (value < 1000000000000) {
            scaled = value / 1000000000;
            suffix = 'B';
        } else {
            scaled = value / 1000000000000;
            suffix = 'T';
        }

        const roundedOne = Math.round(scaled * 10) / 10;

        if (suffix === 'T' && roundedOne < 100) {
            return `${roundedOne.toFixed(1)} ${suffix}`;
        }

        if (roundedOne < 10) {
            return `${roundedOne.toFixed(1)} ${suffix}`;
        }

        return `${Math.round(scaled)} ${suffix}`;
    }

    // Initialize number formatting for stats
    document.querySelectorAll('.stat-number').forEach(element => {
        const targetAttr = element.getAttribute('data-target');
        if (!targetAttr) {
            return;
        }

        const target = Number(targetAttr);
        if (Number.isNaN(target)) {
            element.setAttribute('data-target', 0);
            element.setAttribute('data-formatted', '0');
            return;
        }
        element.setAttribute('data-target', target);
        element.setAttribute('data-formatted', formatNumber(target));
    });

    // Tab switching functionality for species pages
    document.querySelectorAll('.nav-link[data-bs-toggle="pill"]').forEach(tab => {
        tab.addEventListener('shown.bs.tab', function (e) {
            const targetId = e.target.getAttribute('data-bs-target');
            const iframe = document.querySelector(targetId + ' iframe');
            if (iframe && iframe.hasAttribute('data-src')) {
                showLoadingState(targetId);
                setupIframeHandlers(iframe);
                iframe.src = iframe.getAttribute('data-src');
                iframe.removeAttribute('data-src');
            }
        });
    });

    // Preload Dash in the background when the species info page is the default tab.
    const dashboardIframe = document.getElementById('dashboard-iframe');
    if (dashboardIframe && dashboardIframe.hasAttribute('data-src')) {
        window.setTimeout(() => {
            if (!dashboardIframe.hasAttribute('data-src')) {
                return;
            }
            const dashboardPane = dashboardIframe.closest('.tab-pane');
            if (dashboardPane) {
                showLoadingState('#' + dashboardPane.id);
            }
            setupIframeHandlers(dashboardIframe);
            dashboardIframe.src = dashboardIframe.getAttribute('data-src');
            dashboardIframe.removeAttribute('data-src');
        }, 300);
    }

    // Add animation classes on scroll
    const observerOptions = {
        threshold: 0.1,
        rootMargin: '0px 0px -50px 0px'
    };

    const observer = new IntersectionObserver((entries) => {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                entry.target.classList.add('fade-in-up');
            }
        });
    }, observerOptions);

    document.querySelectorAll('.species-card, .feature-card').forEach(card => {
        observer.observe(card);
    });
});

// Enhanced iframe management functions
function initializeIframeManagement() {
    // Setup iframe loading handlers
    document.querySelectorAll('iframe').forEach(iframe => {
        if (!iframe.hasAttribute('data-src')) {
            setupIframeHandlers(iframe);
        }
    });
}

function setupIframeHandlers(iframe) {
    if (iframe.dataset.handlersInitialized === 'true') {
        return;
    }
    iframe.dataset.handlersInitialized = 'true';

    const iframeId = iframe.id;
    
    // Set loading timeout
    const loadingTimeout = setTimeout(() => {
        showErrorState(iframeId);
    }, 15000); // 15 second timeout
    
    iframe.addEventListener('load', function() {
        clearTimeout(loadingTimeout);
        try {
            // Check if iframe content is accessible (same-origin)
            const iframeDoc = iframe.contentDocument || iframe.contentWindow.document;
            showIframeContent(iframeId);
        } catch (e) {
            // Cross-origin content, assume it loaded successfully
            showIframeContent(iframeId);
        }
    });
    
    iframe.addEventListener('error', function() {
        clearTimeout(loadingTimeout);
        showErrorState(iframeId);
    });
}

function showLoadingState(containerId) {
    const container = document.querySelector(containerId);
    if (!container) return;
    
    const loading = container.querySelector('.iframe-loading-state');
    const iframe = container.querySelector('iframe');
    const error = container.querySelector('.iframe-error-state');
    
    if (loading) loading.style.display = 'block';
    if (iframe) iframe.style.display = 'none';
    if (error) error.style.display = 'none';
}

function showIframeContent(iframeId) {
    const iframe = document.getElementById(iframeId);
    if (!iframe) return;
    
    const container = iframe.closest('.iframe-container');
    const loading = container.querySelector('.iframe-loading-state');
    const error = container.querySelector('.iframe-error-state');
    
    if (loading) loading.style.display = 'none';
    if (error) error.style.display = 'none';
    iframe.style.display = 'block';
    
    // Announce to screen readers
    iframe.setAttribute('aria-busy', 'false');
}

function showErrorState(iframeId) {
    const iframe = document.getElementById(iframeId);
    if (!iframe) return;
    
    const container = iframe.closest('.iframe-container');
    const loading = container.querySelector('.iframe-loading-state');
    const error = container.querySelector('.iframe-error-state');
    
    if (loading) loading.style.display = 'none';
    iframe.style.display = 'none';
    if (error) error.style.display = 'block';
    
    // Announce to screen readers
    iframe.setAttribute('aria-busy', 'false');
}

function retryIframe(iframeId) {
    const iframe = document.getElementById(iframeId);
    if (!iframe) return;
    
    const container = iframe.closest('.iframe-container');
    showLoadingState('#' + container.id);
    
    // Force reload
    const currentSrc = iframe.src;
    iframe.src = '';
    setTimeout(() => {
        iframe.src = currentSrc;
        setupIframeHandlers(iframe);
    }, 100);
}

// Add keyboard navigation for tabs
function handleKeyboardNavigation(event) {
    if (event.target.getAttribute('role') === 'tab') {
        const tabs = Array.from(document.querySelectorAll('[role="tab"]'));
        const currentIndex = tabs.indexOf(event.target);
        
        let newIndex;
        if (event.key === 'ArrowRight' || event.key === 'ArrowDown') {
            newIndex = (currentIndex + 1) % tabs.length;
        } else if (event.key === 'ArrowLeft' || event.key === 'ArrowUp') {
            newIndex = (currentIndex - 1 + tabs.length) % tabs.length;
        } else if (event.key === 'Home') {
            newIndex = 0;
        } else if (event.key === 'End') {
            newIndex = tabs.length - 1;
        }
        
        if (newIndex !== undefined) {
            event.preventDefault();
            tabs[newIndex].focus();
            tabs[newIndex].click();
        }
    }
}

document.addEventListener('keydown', handleKeyboardNavigation);

// Search functionality
function initializeSearch() {
    const searchInput = document.getElementById('species-search');
    const clearButton = document.getElementById('clear-search');
    const searchInfo = document.getElementById('search-info');
    const speciesGrid = document.getElementById('species-grid');
    
    if (!searchInput || !speciesGrid) return;
    
    let searchTimeout;
    
    // Debounced search function
    function performSearch() {
        const query = searchInput.value.toLowerCase().trim();
        const allSpecies = speciesGrid.querySelectorAll('[data-species]');
        let visibleCount = 0;
        
        allSpecies.forEach(function(speciesCard) {
            const searchTerms = speciesCard.getAttribute('data-search-terms') || '';
            const isMatch = query === '' || searchTerms.includes(query);
            
            if (isMatch) {
                speciesCard.style.display = 'block';
                speciesCard.classList.remove('search-hidden');
                visibleCount++;
            } else {
                speciesCard.style.display = 'none';
                speciesCard.classList.add('search-hidden');
            }
        });
        
        // Update search info
        if (query === '') {
            searchInfo.textContent = '';
        } else if (visibleCount === 0) {
            searchInfo.textContent = 'No species found matching your search.';
            searchInfo.className = 'search-results-info mt-2 text-warning small';
        } else if (visibleCount === allSpecies.length) {
            searchInfo.textContent = 'Showing all species.';
            searchInfo.className = 'search-results-info mt-2 text-muted small';
        } else {
            searchInfo.textContent = `Showing ${visibleCount} of ${allSpecies.length} species.`;
            searchInfo.className = 'search-results-info mt-2 text-success small';
        }
    }
    
    // Search input handler with debouncing
    searchInput.addEventListener('input', function() {
        clearTimeout(searchTimeout);
        searchTimeout = setTimeout(performSearch, 300);
    });
    
    // Clear search functionality
    if (clearButton) {
        clearButton.addEventListener('click', function() {
            searchInput.value = '';
            performSearch();
            searchInput.focus();
        });
    }
    
    // Show/hide clear button based on input
    searchInput.addEventListener('input', function() {
        if (clearButton) {
            clearButton.style.display = this.value ? 'block' : 'none';
        }
    });
    
    // Initialize clear button state
    if (clearButton) {
        clearButton.style.display = 'none';
    }
    
    // Keyboard shortcuts
    document.addEventListener('keydown', function(e) {
        // Ctrl/Cmd + K to focus search
        if ((e.ctrlKey || e.metaKey) && e.key === 'k') {
            e.preventDefault();
            searchInput.focus();
            searchInput.select();
        }
        
        // Escape to clear search when focused
        if (e.key === 'Escape' && document.activeElement === searchInput) {
            searchInput.value = '';
            performSearch();
            searchInput.blur();
        }
    });
}

// Add CSS for animations and loading states
const style = document.createElement('style');
style.textContent = `
    .spin {
        animation: spin 1s linear infinite;
    }
    
    @keyframes spin {
        from { transform: rotate(0deg); }
        to { transform: rotate(360deg); }
    }
    
    .fade-in-up {
        animation: fadeInUp 0.8s ease forwards;
    }
    
    @keyframes fadeInUp {
        from {
            opacity: 0;
            transform: translateY(30px);
        }
        to {
            opacity: 1;
            transform: translateY(0);
        }
    }
    
    .iframe-loading-state,
    .iframe-error-state {
        background: #f8f9fa;
        border-radius: 15px;
        min-height: 400px;
        display: flex;
        align-items: center;
        justify-content: center;
    }
    
    .spinner-border {
        width: 3rem;
        height: 3rem;
    }
`;
document.head.appendChild(style);

// SequenceServer URL fixing functionality
function initializeSequenceServerFix() {
    // Disabled: direct access to /miniodp/sequenceserver/ works on the running server,
    // so additional URL rewriting inside the iframe is more likely to break
    // SequenceServer's own routing than to help it.
}

function fixSequenceServerUrls(iframe) {
    return;
}

function scheduleUrlCheck(iframe) {
    return;
}
